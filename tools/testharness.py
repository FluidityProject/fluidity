#!/usr/bin/env python3

import json
import re
import subprocess
import sys
from argparse import ArgumentParser
from os import environ, sched_getaffinity
from pathlib import Path
from textwrap import dedent, indent
from time import monotonic, sleep
from xml.etree.ElementTree import parse

try:
    from importlib.metadata import version

    from junit_xml import TestCase, TestSuite

    assert (
        float(version("junit_xml")) >= 1.9
    ), "ERROR: junit_xml version must be at least 1.9 — please update."
except ImportError:  # Provide dummy classes in the absence of junit_xml

    class TestSuite(object):
        def __init__(self, *args, **kwargs):
            self.test_cases = []

        def to_file(self, *args, **kwargs):
            print("WARNING: junit_xml required to produce an XML output file.")

    class TestCase(object):
        def __init__(self, *args, **kwargs):
            pass

        def add_error_info(self, *args, **kwargs):
            pass

        def add_failure_info(self, *args, **kwargs):
            pass

        def add_skipped_info(self, *args, **kwargs):
            pass


def filter_tests(xml_files, test_suite):
    """Iterate through all found XML files and retain only the ones that match the
    program input arguments. Update the JUnit record with skipped tests."""
    tests = {}
    for xml_file in xml_files:
        # Obtain basic information about the test
        parsed_xml = parse(xml_file).getroot()
        if parsed_xml.tag != "testproblem":
            print(f"Skipping {parsed_xml} — root tag differs from 'testproblem'.")
            continue
        prob_def = parsed_xml.find("problem_definition")
        prob_length = prob_def.get("length")
        prob_nprocs = int(prob_def.get("nprocs"))
        # Manually insert the required amount of cores for parallel tests
        prob_command = prob_def.find("command_line").text.replace(
            "mpiexec", f"mpiexec -n {prob_nprocs}"
        )
        # Indicate if the test belongs to 'examples' or 'tests' and specify its length.
        prob_classname = f"{xml_file.parts[-3]}.{prob_length}"
        try:
            xml_tags = parsed_xml.find("tags").text.split()
        except AttributeError:  # If no tags are present, assign an empty list
            xml_tags = []
        # Conditions to exclude a test
        length_condition = (
            False if args.length is None else prob_length not in args.length
        )
        nprocs_condition = (args.exec_type == "parallel" and prob_nprocs <= 1) or (
            args.exec_type == "serial" and prob_nprocs != 1
        )
        omitted_tags_condition = set(xml_tags).intersection(args.omit_tags)
        required_tags_condition = set(args.tags).difference(xml_tags)
        # Skip the test if any of the conditions is not met and update the JUnit record
        if (
            length_condition
            or nprocs_condition
            or omitted_tags_condition
            or required_tags_condition
        ):
            test_case = TestCase(
                name=xml_file.stem, classname=prob_classname, status="Skipped"
            )
            test_case.add_skipped_info(message="Test suite conditions unmet.")
            test_suite.test_cases.append(test_case)
            continue
        # Populate dictionary with test information
        tests[xml_file] = {
            "id": prob_classname,
            "length": prob_length,
            "nprocs": prob_nprocs,
            "command": prob_command,
            "variables": {},
            "pass_tests": {},
            "warn_tests": {},
            "elapsed_time": 0,
            "stdout": "\n",
            "stderr": "\n",
        }
        for var in parsed_xml.iter("variable"):
            assert var.get("language").strip() == "python"
            tests[xml_file]["variables"][var.get("name")] = dedent(var.text)
        for test in parsed_xml.find("pass_tests"):
            assert test.get("language").strip() == "python"
            tests[xml_file]["pass_tests"][test.get("name")] = dedent(test.text)
        try:
            for test in parsed_xml.find("warn_tests"):
                assert test.get("language").strip() == "python"
                tests[xml_file]["warn_tests"][test.get("name")] = dedent(test.text)
        except TypeError:  # This test does not contain Warn Tests
            pass
    return tests


def gather_tests():
    """Look for tests given the program input arguments."""
    if args.file:  # Specific test requested
        xml_files = [Path(args.file)]
    elif args.from_file:  # Specific list of tests requested
        with open(args.from_file, "r") as fid:
            xml_files = [Path(test_name.rstrip()) for test_name in fid]
    else:  # Gather all XML files that can be found
        test_paths = ["examples", "tests", "longtests"]
        xml_files = [
            xml_file
            for test_path in test_paths
            if (fluidity_root / test_path).exists()
            for xml_file in (fluidity_root / test_path).glob("*/*.xml")
        ]
    return xml_files


def generate_python_test_string(test_type, assertion_status, test, test_str):
    """Encapsulate within a try/except structure the Python code to execute to properly
    catch potential exceptions and obtain tracebacks."""
    return f"""
print("# {test_type} Test: {test}")
try:
{indent(test_str.strip(), "    ")}
    print("--- Success")
except AssertionError:
    from sys import stderr
    from traceback import print_exc
    print("--- {assertion_status} !!!")
    print_exc()
    print("~~~ End of Traceback ~~~", file=stderr)
except Exception:
    from sys import stderr
    from traceback import print_exc
    print("--- Error !!!")
    print_exc()
    print("~~~ End of Traceback ~~~", file=stderr)"""


def poll_processes(
    running_procs,
    core_counter,
    serial,
    return_list,
    task_string,
    process_interpreter,
    test_suite,
):
    """Check if running processes have terminated and deal with results."""
    proc_status = [
        tests[running_xml]["running_proc"].poll() for running_xml in running_procs
    ]
    # Reversed to be able to remove items from lists without skipping any item
    for status, running_xml in zip(reversed(proc_status), reversed(running_procs)):
        if status is not None:  # Test has terminated
            current_test = tests[running_xml]
            # Measure an upper bound for the test elapsed time
            current_test["elapsed_time"] += monotonic() - current_test["create_time"]
            # Update objects that keep track of the current activity
            task_ncores = 1 if serial else current_test["nprocs"]
            core_counter -= task_ncores
            running_procs.remove(running_xml)
            # Recover standard ouput and error streams
            stdout = "".join([line for line in read_stream(running_xml, "stdout_file")])
            stderr = "".join([line for line in read_stream(running_xml, "stderr_file")])
            current_test["stdout"] += f"{task_string}\n{stdout}\n"
            current_test["stderr"] += f"{task_string}\n{stderr}\n"
            # Deal with successful processes
            if "No rule to make target 'input'" in stderr or (
                status == 0
                and all(
                    x not in stderr for x in ["Traceback", "/bin/bash", "MPI_ABORT"]
                )
            ):
                return_list.append(running_xml)
            else:  # Deal with errors
                process_error(
                    running_xml, process_interpreter, stdout, stderr, test_suite
                )
    sleep(0.1)  # Avoid high-intensity polling
    return core_counter


def process_error(test_xml, process_interpreter, stdout, stderr, test_suite):
    """Process tests that did not complete successfully and update the JUnit record
    accordingly."""
    # Add an entry to the XML parser
    test_case = TestCase(
        name=test_xml.stem,
        classname=tests[test_xml]["id"],
        elapsed_sec=tests[test_xml]["elapsed_time"],
        stdout=tests[test_xml]["stdout"],
        stderr=tests[test_xml]["stderr"],
        allow_multiple_subelements=True,
    )
    # Errors generated in make input or while running command(s)
    if process_interpreter == "Shell":
        error_list.append(test_xml.stem)
        test_case.status = "Error"
        test_case.add_error_info(
            message=f"{test_xml.stem}: Shell script failed.", output=f"\n{stderr}"
        )
        if tests[test_xml]["running_proc"].returncode == 0:
            tests[test_xml]["running_proc"].returncode = "Failed command"
    # Errors/Failures generated in Python variables/tests
    elif process_interpreter == "Python":
        python_lines = tests[test_xml]["running_proc"].args[-1].split("\n")
        # Error while parsing Python test string
        if any(
            error in stderr for error in ["IndentationError", "SyntaxError", "TabError"]
        ):
            error_list.append(test_xml.stem)
            test_case.status = "Error"
            test_case.add_error_info(
                message=f"{test_xml.stem}: Parsing failed", output=f"\n{stderr}"
            )
        # Failure(s) within actual test(s)
        elif any(f"# {kind} Test" in stdout for kind in ["Pass", "Warn"]):
            # Identify unsuccessful Python tests through stdout
            regex = (
                r"(?<=^#\s)"
                ".+"
                "(?=(\n^[^#].+|\n^#(?! Pass Test| Warn Test).+){0,}"
                "\n--- (Error|Failure|Warning))"
            )
            failed_tests = [
                match.group().split(": ", maxsplit=1)
                for match in re.finditer(regex, stdout, re.MULTILINE)
            ]
            # Split stderr into individual tracebacks
            regex = "(^Traceback).+(\n.+)+?(?=\n~~~ End of Traceback ~~~)"
            tracebacks = [
                match.group() for match in re.finditer(regex, stderr, re.MULTILINE)
            ]
            flag_pass, flag_warn = True, True
            for traceback, (test_type, test_name) in zip(tracebacks, failed_tests):
                if flag_pass and test_type == "Pass Test":
                    failure_list.append(test_xml.stem)
                    test_case.status = "Failure"
                    flag_pass = False
                elif flag_warn and test_type == "Warn Test":
                    warning_list.append(test_xml.stem)
                    test_case.status = "Warning" if flag_pass else "Failure and Warning"
                    flag_warn = False
                line_nb = re.search(
                    '(?<="<string>", line )[0-9]+(?=, in <module>)', traceback
                ).group()
                python_error_line = python_lines[int(line_nb) - 1][4:]
                traceback = traceback.replace(
                    f'File "<string>", line {line_nb}, in <module>',
                    f"\n  Caught exception at '{python_error_line.strip()}'\n",
                )
                failure_type = "failure" if "Pass" in test_type else "warning"
                test_case.add_failure_info(
                    message=f"{test_xml.stem}: Test '{test_name}' failed",
                    output=f"\n{traceback}",
                    failure_type=failure_type,
                )
        else:  # Error within variable
            error_list.append(test_xml.stem)
            test_case.status = "Error"
            try:
                line_nb = re.search(
                    '(?<="<string>", line )[0-9]+(?=, in <module>)', stderr
                ).group()
                python_error_line = python_lines[int(line_nb) - 1]
                traceback = stderr.replace(
                    f'File "<string>", line {line_nb}, in <module>',
                    f"\n  Caught exception at '{python_error_line.strip()}'\n",
                )
            except AttributeError:
                traceback = stderr
            var_name = list(
                re.finditer(r"(?<=^#\sVariable:\s).+", stdout, re.MULTILINE)
            )[-1].group()
            test_case.add_error_info(
                message=f"{test_xml.stem}: Variable '{var_name}' failed.",
                output=f"\n{traceback}",
            )
        if tests[test_xml]["running_proc"].returncode == 0:
            tests[test_xml]["running_proc"].returncode = "Python exception"
    test_suite.test_cases.append(test_case)
    # Print relevant information regarding the error
    print(f"\n* ERROR: {test_xml.stem} exited with a non-zero exit code.")
    print(f"* Exit status: {tests[test_xml]['running_proc'].returncode}")
    print(f"* Output:\n{stdout.strip()}")
    print(f"* Stderr output:\n{stderr}")


def read_stream(test_xml, stream_key):
    """Read content of provided stream."""
    # Set file object’s position to the beginning of the file
    tests[test_xml][stream_key].seek(0)
    stream = tests[test_xml][stream_key].readlines()
    if bool(stream) is False:  # Check if stream is empty
        stream = ["No output generated.\n"]
    tests[test_xml][stream_key].close()  # Close file object
    return stream


def run_tasks(
    task_string, tests_list, serial, process_interpreter, task_function, test_suite
):
    """Iterate through the test list and execute the provided function in parallel."""
    return_list, running_procs, core_counter = [], [], 0
    print(task_string)
    # Start tasks whilst there are queuing tests
    while tests_list:
        # Remove and obtain a test from the end of the queue
        test_xml = tests_list.pop()
        current_test = tests[test_xml]
        # Determine the amount of cores required
        task_ncores = 1 if serial else current_test["nprocs"]
        # Check if there are sufficient resources available; allow tests requesting a
        # higher number of cores than available to run by oversubscribing nodes (MPI),
        # with only one oversubscribed test at a time
        if core_counter >= core_avail or (
            core_counter + task_ncores > core_avail and task_ncores <= core_avail
        ):
            # Check if some processes have terminated
            core_counter = poll_processes(
                running_procs,
                core_counter,
                serial,
                return_list,
                task_string,
                process_interpreter,
                test_suite,
            )
            # Check if sufficient resources have become available
            if core_counter >= core_avail or (
                core_counter + task_ncores > core_avail and task_ncores <= core_avail
            ):
                # Re-insert the test at the beginning of the queue
                tests_list.insert(0, test_xml)
                # Skip the remainder of the loop iteration
                continue
        print(f"\t-> New test: {test_xml.stem}")
        # Open streams to re-direct stdout and stderr
        current_test["stdout_file"] = open(test_xml.parent / "stdout", "w+")
        current_test["stderr_file"] = open(test_xml.parent / "stderr", "w+")
        # Register the starting time of the task
        current_test["create_time"] = monotonic()
        # Submit task and update objects that keep track of the current load
        task_function(test_xml)
        core_counter += task_ncores
        running_procs.append(test_xml)
    # Once the queue is empty, wait for the processes that are still running
    while running_procs:
        core_counter = poll_processes(
            running_procs,
            core_counter,
            serial,
            return_list,
            task_string,
            process_interpreter,
            test_suite,
        )

    # Check that the objects which keep track of the current load have come back to
    # their nominal value
    assert core_counter == 0 and bool(running_procs) is False

    return return_list


def set_environment_variable(env_var, env_path):
    """Set or prepend to the requested environment variable."""
    try:
        environ[env_var] = f"{env_path}:{environ[env_var]}"
    except KeyError:  # If the environment variable does not exist, create it
        environ[env_var] = str(env_path)


def task_make_input(test_xml):
    """Execute `make input`."""
    tests[test_xml]["running_proc"] = subprocess.Popen(
        ["make", "input"],
        cwd=test_xml.parent,
        encoding="utf-8",
        stdout=tests[test_xml]["stdout_file"],
        stderr=tests[test_xml]["stderr_file"],
    )


def task_run_commands(test_xml):
    """Execute test instructions."""
    tests[test_xml]["running_proc"] = subprocess.Popen(
        tests[test_xml]["command"],
        cwd=test_xml.parent,
        shell=True,
        encoding="utf-8",
        stdout=tests[test_xml]["stdout_file"],
        stderr=tests[test_xml]["stderr_file"],
    )


def task_run_tests(test_xml):
    """Calculate Python variables specific to the test and assess their values."""
    # Join variable strings together
    var_string = "\n".join(
        [
            f"print('# Variable: {var}')\n{var_str.strip()}"
            for var, var_str in tests[test_xml]["variables"].items()
        ]
    )
    # Join pass-test strings together
    pass_string = "\n".join(
        [
            generate_python_test_string("Pass", "Failure", test, test_str)
            for test, test_str in tests[test_xml]["pass_tests"].items()
        ]
    )
    # Join warn-test strings together
    warn_string = "\n".join(
        [
            generate_python_test_string("Warn", "Warning", test, test_str)
            for test, test_str in tests[test_xml]["warn_tests"].items()
        ]
    )
    # Join all three strings
    test_string = "\n".join([var_string, pass_string, warn_string])
    tests[test_xml]["running_proc"] = subprocess.Popen(
        [sys.executable, "-c", test_string],
        cwd=test_xml.parent,
        encoding="utf-8",
        stdout=tests[test_xml]["stdout_file"],
        stderr=tests[test_xml]["stderr_file"],
    )


parser = ArgumentParser(description="Fluidity Test Harness")
parser.add_argument(
    "--clean", action="store_true", help="call `make clean` for each test found"
)
parser.add_argument(
    "-e",
    "--exec-type",
    default="any",
    choices=["serial", "parallel"],
    help="specify which kind of tests to run; choose either serial or parallel",
    metavar="TYPE",
)
parser.add_argument("-f", "--file", help="run a single test — expects an XML filepath")
parser.add_argument(
    "--from-file",
    help="path to a file where to read which tests to run — one XML filepath per line",
    metavar="FILE",
)
parser.add_argument(
    "--just-list",
    nargs="?",
    const=True,
    default=False,
    help="print which tests were found and save the list to a JSON file if provided",
    metavar="FILE",
)
parser.add_argument(
    "--just-test",
    action="store_true",
    help="execute Python instructions without re-running test core commands",
)
parser.add_argument(
    "-l",
    "--length",
    action="append",
    choices=["vshort", "short", "medium", "long", "vlong"],
    help="test length(s) to be run; choose from vshort, short, medium, long or vlong",
)
parser.add_argument(
    "-n",
    "--ncores",
    type=int,
    help="number of logical cores to target",
    metavar="CORES",
)
parser.add_argument(
    "-o",
    "--omit-tags",
    action="extend",
    nargs="*",
    default=[],
    help="tags identifying which tests to exclude",
    metavar="TAG",
)
parser.add_argument(
    "-t",
    "--tags",
    action="extend",
    nargs="*",
    default=[],
    help="tags identifying which tests to run",
    metavar="TAG",
)
parser.add_argument("-v", "--valgrind", action="store_true", help="enable Valgrind")
parser.add_argument("-x", "--xml-output", help="XML output filename", metavar="OUTPUT")
args = parser.parse_args()

# Obtain path to the root of fluidity's directory
fluidity_root = Path(sys.argv[0]).resolve().parent.parent

set_environment_variable("PATH", fluidity_root / "bin")
set_environment_variable("PATH", fluidity_root / "libspud" / "bin")
set_environment_variable("PYTHONPATH", fluidity_root / "python")
set_environment_variable("PYTHONPATH", fluidity_root / "libspud" / "diamond")
set_environment_variable("LD_LIBRARY_PATH", fluidity_root / "lib")

print(
    f"""*** Test criteria
\t-> length: {"any" if args.length is None else " ".join(args.length)}
\t-> parallel: {args.exec_type}
\t-> tags to include: {args.tags}
\t-> tags to exclude: {args.omit_tags}
"""
)

xml_files = gather_tests()
assert xml_files, "No tests were found."
test_suite = TestSuite("Test Harness")
tests = filter_tests(xml_files, test_suite)
assert tests, "No tests matched the provided test criteria."

if args.just_list:
    print("*** Found tests that match the input criteria")
    for test_xml in sorted(tests.keys()):
        print(f"\t-> {test_xml.stem}")
    nb_tests = len(tests.keys())
    print(f"{nb_tests} test{'s' if nb_tests > 1 else ''} found.")
    if isinstance(args.just_list, str):
        with open(args.just_list, "w") as fid:
            json.dump([test.name for test in tests.keys()], fid)
    raise SystemExit
elif args.clean:
    print("*** Cleaning")
    for test_xml in tests.keys():
        print(f"\t-> {test_xml.stem}: Calling `make clean`.")
        try:
            subprocess.run(
                ["make", "clean"],
                cwd=test_xml.parent,
                encoding="utf-8",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
        except subprocess.CalledProcessError as test_error:
            if "No rule to make target 'clean'" in test_error.stderr:
                pass
            else:
                raise test_error
    print("Cleaning done.")
    raise SystemExit

fluidity_version = subprocess.run(
    ["fluidity", "-V"],
    check=True,
    encoding="utf-8",
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
)
print(
    f"""{"-" * 80}
Output of "fluidity -V"
{fluidity_version.stderr.strip()}
{"-" * 80}
"""
)
if args.valgrind:
    # Make sure Fluidity has been compiled in debug mode
    assert (
        "debugging" in fluidity_version.stderr
    ), "Please compile Fluidity in debug mode to use valgrind."
    print(
        f"""{"-" * 80}
I see you are using Valgrind!
Keep the following in mind:
- The log file will be produced in the directory containing the tests.
- Valgrind typically takes O(100) times as long. I hope your test is short.
{"-" * 80}
"""
    )
    # Prepend valgrind command
    for test_xml in tests.keys():
        # What about when command contains commands delimited by ";"?
        tests[test_xml][
            "command"
        ] = f"""valgrind --tool=memcheck \
--leak-check=full -v --show-reachable=yes --num-callers=8 --error-limit=no \
--log-file=test.log {tests[test_xml]["command"]}"""

core_avail = len(sched_getaffinity(0))
if args.ncores is not None:
    core_avail = min(args.ncores, core_avail)
error_list, failure_list, warning_list = [], [], []

tests_list = list(tests.keys())
if not args.just_test:
    tests_list = run_tasks(
        "*** Executing 'make input'",
        tests_list,
        True,
        "Shell",
        task_make_input,
        test_suite,
    )
    if tests_list:
        tests_list = run_tasks(
            "*** Executing test commands from XML files",
            tests_list,
            False,
            "Shell",
            task_run_commands,
            test_suite,
        )
if tests_list:
    tests_list = run_tasks(
        "*** Executing Python tests",
        tests_list,
        True,
        "Python",
        task_run_tests,
        test_suite,
    )

for test_xml in tests_list:
    test_case = TestCase(
        name=test_xml.stem,
        classname=tests[test_xml]["id"],
        elapsed_sec=tests[test_xml]["elapsed_time"],
        status="Success",
    )
    test_case.stdout = tests[test_xml]["stdout"]
    test_case.stderr = tests[test_xml]["stderr"]
    test_suite.test_cases.append(test_case)

if args.xml_output:
    with open(args.xml_output, "w") as fid:
        test_suite.to_file(fid, [test_suite])

if any([error_list, failure_list, warning_list]):
    if error_list:
        print("Summary of test problems that produced errors:")
        for test in error_list:
            print(f"-> {test}")
    if failure_list:
        print("Summary of test problems that produced failures:")
        for test in failure_list:
            print(f"-> {test}")
    if warning_list:
        print("Summary of test problems that produced warnings:")
        for test in warning_list:
            print(f"-> {test}")
else:
    print("Test suite completed successfully.")
