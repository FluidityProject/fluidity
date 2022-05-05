#!/usr/bin/env python3

import json
import re
import subprocess
import sys
from argparse import ArgumentParser
from itertools import chain
from os import environ, sched_getaffinity
from pathlib import Path
from textwrap import dedent, indent
from time import monotonic, sleep
from xml.etree.ElementTree import parse

try:
    from junit_xml import TestCase, TestSuite

    if sys.version_info.major == 3 and sys.version_info.minor >= 8:
        from importlib.metadata import version

        assert (
            float(version("junit_xml")) >= 1.9
        ), "ERROR: junit_xml version must be at least 1.9 - please update."
    else:  # To be dropped when Python >= 3.8 becomes mainstream
        from pkg_resources import get_distribution

        assert (
            float(get_distribution("junit_xml").version) >= 1.9
        ), "ERROR: junit_xml version must be at least 1.9 - please update."
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


def filter_tests(xml_files):
    tests = {}
    # Iterate through all found XML files and keep only the ones that match the
    # program input arguments
    for xml_file in xml_files:
        # Obtain basic information about the test
        parsed_xml = parse(xml_file).getroot()
        assert parsed_xml.tag == "testproblem", xml_file
        prob_def = parsed_xml.find("problem_definition")
        prob_length = prob_def.get("length")
        prob_nprocs = int(prob_def.get("nprocs"))
        # Manually insert the required amount of cores for parallel tests
        prob_command = prob_def.find("command_line").text.replace(
            "mpiexec", f"mpiexec -n {prob_nprocs}"
        )
        prob_classname = f"{xml_file.parts[-3]}.{prob_length}"
        try:
            xml_tags = parsed_xml.find("tags").text.split()
        except AttributeError:  # If no tags are present, assign an empty list
            xml_tags = []
        # Define conditions the test must meet
        # Do we want to keep this condition, or do we want to have "any"
        # correspond to all tests when length is not given (and only run the
        # provided lengths when length is given)?
        length_condition = (
            args.length == "any" and prob_length in ["long", "vlong"]
        ) or (args.length != "any" and prob_length not in args.length)
        nprocs_condition = (args.parallel == "parallel" and prob_nprocs <= 1) or (
            args.parallel == "serial" and prob_nprocs != 1
        )
        excluded_tags_condition = set(xml_tags).intersection(args.exclude_tags)
        required_tags_condition = set(args.tags).difference(xml_tags)
        # Skip the test if any of the conditions is not met
        if (
            length_condition
            or nprocs_condition
            or excluded_tags_condition
            or required_tags_condition
        ):
            xml_entry = TestCase(
                name=xml_file.stem, classname=prob_classname, status="Skipped"
            )
            xml_entry.add_skipped_info(message="Test suite conditions unmet.")
            xml_parser.test_cases.append(xml_entry)
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
    if args.file:  # Specific test requested
        xml_files = [fluidity_source.rglob(args.file)]
    elif args.from_file:  # Specific list of tests requested
        with open(args.from_file, "r") as fid:
            xml_files = [fluidity_source.rglob(test_name.rstrip()) for test_name in fid]
    else:  # Gather all XML files that can be found
        test_paths = ["examples", "tests", "longtests"]
        xml_files = [
            list((fluidity_source / test_path).rglob("*/*.xml"))
            for test_path in test_paths
            if (fluidity_source / test_path).exists()
        ]
    return list(chain(*xml_files))


def generate_string(test_type, assertion_status, test, test_str):
    # Encapsulate Python code to execute within a try/except structure to
    # properly catch potential exceptions and obtain tracebacks
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
    oversubscribe,
):
    # Check if running processes have terminated
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
            task_nprocs = 1 if serial else current_test["nprocs"]
            core_counter -= task_nprocs
            if task_nprocs > core_avail:
                oversubscribe = False
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
                process_error(running_xml, process_interpreter, stdout, stderr)
    sleep(0.1)  # Avoid high-intensity polling
    return core_counter, oversubscribe


def process_error(test_xml, process_interpreter, stdout, stderr):
    # Add an entry to the XML parser
    xml_entry = TestCase(
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
        xml_entry.status = "Error"
        xml_entry.add_error_info(
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
            xml_entry.status = "Error"
            xml_entry.add_error_info(
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
                    xml_entry.status = "Failure"
                    flag_pass = False
                elif flag_warn and test_type == "Warn Test":
                    warning_list.append(test_xml.stem)
                    xml_entry.status = "Warning" if flag_pass else "Failure and Warning"
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
                xml_entry.add_failure_info(
                    message=f"{test_xml.stem}: Test '{test_name}' failed",
                    output=f"\n{traceback}",
                    failure_type=failure_type,
                )
        else:  # Error within variable
            error_list.append(test_xml.stem)
            xml_entry.status = "Error"
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
            xml_entry.add_error_info(
                message=f"{test_xml.stem}: Variable '{var_name}' failed.",
                output=f"\n{traceback}",
            )
        if tests[test_xml]["running_proc"].returncode == 0:
            tests[test_xml]["running_proc"].returncode = "Python exception"
    xml_parser.test_cases.append(xml_entry)
    # Print relevant information regarding the error
    print(f"\n* ERROR: {test_xml.stem} exited with a non-zero exit code.")
    print(f"* Exit status: {tests[test_xml]['running_proc'].returncode}")
    print(f"* Output:\n{stdout.strip()}")
    print(f"* Stderr output:\n{stderr}")


def read_stream(test_xml, stream_key):
    # Set file object’s position to the beginning of the file
    tests[test_xml][stream_key].seek(0)
    stream = tests[test_xml][stream_key].readlines()
    if bool(stream) is False:  # Check if stream is empty
        stream = ["No output generated.\n"]
    tests[test_xml][stream_key].close()  # Close file object
    return stream


def run_tasks(task_string, tests_list, serial, process_interpreter, task_function):
    return_list, running_procs, core_counter = [], [], 0
    print(task_string)
    oversubscribe = False
    # Start tasks whilst there are queuing tests
    while tests_list:
        # Remove and obtain a test from the end of the queue
        test_xml = tests_list.pop()
        current_test = tests[test_xml]
        # Determine the amount of cores required
        task_nprocs = 1 if serial else current_test["nprocs"]
        # Check if there are sufficient resources available
        if core_counter + task_nprocs > core_avail and (
            task_nprocs <= core_avail or oversubscribe
        ):
            # Check if some processes have terminated
            core_counter, oversubscribe = poll_processes(
                running_procs,
                core_counter,
                serial,
                return_list,
                task_string,
                process_interpreter,
                oversubscribe,
            )
            # Check if sufficient resources have become available
            if core_counter + task_nprocs > core_avail and (
                task_nprocs <= core_avail or oversubscribe
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
        core_counter += task_nprocs
        running_procs.append(test_xml)
        # Allow tests requesting a higher number of cores than available to
        # run by oversubscribing nodes (MPI)
        if task_nprocs > core_avail:
            oversubscribe = True
    # Once the queue is empty, wait for the processes that are still running
    while running_procs:
        core_counter, oversubscribe = poll_processes(
            running_procs,
            core_counter,
            serial,
            return_list,
            task_string,
            process_interpreter,
            oversubscribe,
        )

    # Check that the objects which keep track of the current load have come
    # back to their nominal value
    assert core_counter == 0 and bool(running_procs) is False
    assert oversubscribe is False

    return return_list


def set_environment_variable(env_var, env_path):
    try:  # Check if the the environment variable already contains the path
        if str(env_path) not in environ[env_var]:
            environ[env_var] = f"{env_path}:" + environ[env_var]
    except KeyError:  # If the environment variable does not exist, create it
        environ[env_var] = str(env_path)


def task_make_input(test_xml):
    tests[test_xml]["running_proc"] = subprocess.Popen(
        ["make", "input"],
        cwd=test_xml.parent,
        encoding="utf-8",
        stdout=tests[test_xml]["stdout_file"],
        stderr=tests[test_xml]["stderr_file"],
    )


def task_run_commands(test_xml):
    tests[test_xml]["running_proc"] = subprocess.Popen(
        tests[test_xml]["command"],
        cwd=test_xml.parent,
        shell=True,
        encoding="utf-8",
        stdout=tests[test_xml]["stdout_file"],
        stderr=tests[test_xml]["stderr_file"],
    )


def task_run_tests(test_xml):
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
            generate_string("Pass", "Failure", test, test_str)
            for test, test_str in tests[test_xml]["pass_tests"].items()
        ]
    )
    # Join warn-test strings together
    warn_string = "\n".join(
        [
            generate_string("Warn", "Warning", test, test_str)
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
    "-c",
    "--clean",
    action="store_true",
    help="call 'make clean' for each included test",
)
# action="extend", nargs=1 can be used with Python >= 3.8
parser.add_argument(
    "-e",
    "--exclude-tags",
    metavar="",
    action="append",
    default=[],
    help="tags indicating tests that are excluded",
)
parser.add_argument(
    "-f", "--file", metavar="", help="single test to run - expects XML filename"
)
# action="extend", nargs=1 can be used with Python >= 3.8
parser.add_argument(
    "-l",
    "--length",
    metavar="",
    action="append",
    help="""test length to be run; must be either vshort, short, medium, long or
vlong""",
)
parser.add_argument(
    "-n", "--nprocs", type=int, metavar="", help="targeted number of cores"
)
parser.add_argument(
    "-p",
    "--parallelism",
    dest="parallel",
    default="any",
    metavar="",
    help="""parallelism of problem; must be either serial, parallel or any - defaults
to '%(default)s'""",
)
parser.add_argument(
    "-t",
    "--tags",
    metavar="",
    action="extend",
    nargs=1,
    default=[],
    help="tags indicating tests that are included",
)
parser.add_argument("-v", "--valgrind", action="store_true", help="enable Valgrind")
parser.add_argument("-x", "--xml-output", metavar="", help="XML output filename")
parser.add_argument(
    "--from-file",
    metavar="",
    help="""path to a file containing a list of tests to run (one *.xml filename per
line)""",
)
parser.add_argument(
    "--github",
    metavar="",
    help="""filename (expects JSON) to store the tests to run on GitHub Actions""",
)
parser.add_argument(
    "--just-list",
    action="store_true",
    help="""prints the list of tests that would be included in the run given the
provided arguments""",
)
args = parser.parse_args()

if args.length is None:
    args.length = "any"
elif set(args.length).difference(["vshort", "short", "medium", "long", "vlong"]):
    parser.error(
        """Specify length as either of vshort, short, medium, long, or vlong."""
    )
if args.parallel not in ["serial", "parallel", "any"]:
    parser.error("Specify parallelism as either of serial, parallel or any.")

# Obtain path to the root of fluidity's directory
fluidity_source = Path("@CMAKE_SOURCE_DIR@")
fluidity_build = Path("@CMAKE_BINARY_DIR@")

set_environment_variable("OMPI_MCA_rmaps_base_oversubscribe", 1)
set_environment_variable("PATH", fluidity_build / "bin")

print(
    f"""*** Test criteria
\t-> length: {args.length}
\t-> parallel: {args.parallel}
\t-> tags to include: {args.tags}
\t-> tags to exclude: {args.exclude_tags}
"""
)

xml_files = gather_tests()
assert xml_files, "No tests were found."
xml_parser = TestSuite("Test Harness")
tests = filter_tests(xml_files)
assert tests, "No tests matched the provided test criteria."

if args.just_list:
    print("*** Found tests that match the input criteria")
    for test_xml in sorted(tests.keys()):
        print(f"\t-> {test_xml.stem}")
    print(f"{len(tests.keys())} tests found.")
    if args.github:
        with open(args.github, "w") as fid:
            json.dump([test.name for test in tests.keys()], fid)
elif args.clean:
    print("*** Cleaning")
    for test_xml in tests.keys():
        print(f"\t-> {test_xml.stem}: Calling 'make clean'.")
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
else:
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
    if args.nprocs is not None:
        core_avail = min(args.nprocs, core_avail)
    error_list, failure_list, warning_list = [], [], []

    tests_list = run_tasks(
        "*** Executing 'make input'", list(tests.keys()), True, "Shell", task_make_input
    )
    if tests_list:
        tests_list = run_tasks(
            "*** Executing test commands from XML files",
            tests_list,
            False,
            "Shell",
            task_run_commands,
        )
    if tests_list:
        tests_list = run_tasks(
            "*** Executing Python tests", tests_list, True, "Python", task_run_tests
        )

    for test_xml in tests_list:
        xml_entry = TestCase(
            name=test_xml.stem,
            classname=tests[test_xml]["id"],
            elapsed_sec=tests[test_xml]["elapsed_time"],
            status="Success",
        )
        xml_entry.stdout = tests[test_xml]["stdout"]
        xml_entry.stderr = tests[test_xml]["stderr"]
        xml_parser.test_cases.append(xml_entry)

    if args.xml_output:
        with open(args.xml_output, "w") as fid:
            xml_parser.to_file(fid, [xml_parser])

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
