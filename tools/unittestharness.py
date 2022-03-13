#!/usr/bin/env python3

# Pass, Warn and Fail are defined in femtools/Unittest_tools.F90
# Make sure any changes there are reflected here

import argparse
from collections import Counter
from inspect import cleandoc
import os
from pathlib import Path
import re
from shutil import rmtree
import subprocess
import sys
from time import process_time


def display_results(tests_results, error_list, skip_list):
    print(cleandoc(f"""
                   ###########
                   # Results #
                   ###########

                   \tPasses: {sum(tests_results['Pass'].values())}\n
                   \tWarnings: {sum(tests_results['Warn'].values())}
                   \t\t{list(tests_results['Warn'].keys())}\n
                   \tFailures: {sum(tests_results['Fail'].values())}
                   \t\t{list(tests_results['Fail'].keys())}\n
                   \tErrors: {len(error_list)}
                   \t\t{error_list}\n
                   \tSkipped: {len(skip_list)}
                   \t\t{skip_list}"""))


def unittest_harness_no_output(tests):
    tests_results = {'Pass': Counter(), 'Warn': Counter(), 'Fail': Counter()}
    error_list, skip_list = [], []

    for test in tests:
        print(f'\t-> New test: {test.name}')

        if (tests_dir / test.name).is_file() is False:
            print(f'WARNING: {test.name} not found')
            skip_list.append(test.name)
            continue

        try:
            test_proc = subprocess.run(
                './' + test.name, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=True, cwd=test.parent,
                check=True, encoding='utf-8')
        except subprocess.CalledProcessError as test_error:
            print(f'ERROR: {test.name} exited with a non-zero exit code.')
            print(f'Exit status: {test_error.returncode}')
            print(f'Output: {test_error.output}')
            print(f'Stderr output: {test_error.stderr}')
            error_list.append(test.name)
            continue

        for test_output in test_proc.stdout.splitlines():
            try:
                tests_results[test_output[:4]][test.name] += 1
            except KeyError:
                print(f'\t\t\t{test_output.lstrip()}')
                continue

            print(f'\t\t{test_output}')

    display_results(tests_results, error_list, skip_list)


def unittest_harness(tests, xml_outfile):
    xml_parser = TestSuite('unittest_harness')
    tests_results = {'Pass': Counter(), 'Warn': Counter(), 'Fail': Counter()}
    error_list, skip_list = [], []

    for test in tests:
        print(f'\t-> New test: {test.name}')

        if (tests_dir / test.name).is_file() is False:
            print(f'WARNING: {test.name} not found')
            skip_list.append(test.name)
            xml_entry = TestCase('', classname=test.name,
                                 elapsed_sec='None', status='Skipped')
            xml_entry.add_skipped_info(message='Not found')
            xml_parser.test_cases.append(xml_entry)
            continue

        test_id = f'unittest.{test.resolve().parts[-3]}.{test.name}'

        try:
            begin = process_time()
            test_proc = subprocess.run(
                './' + test.name, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=True, cwd=test.parent,
                check=True, encoding='utf-8')
            wall_time = process_time() - begin
        except subprocess.CalledProcessError as test_error:
            print(f'ERROR: {test.name} exited with a non-zero exit code.')
            print(f'Exit status: {test_error.returncode}')
            print(f'Output: {test_error.output}')
            print(f'Stderr output: {test_error.stderr}')
            error_list.append(test.name)
            xml_entry = TestCase('', classname=test_id,
                                 elapsed_sec=process_time() - begin,
                                 status='Error')
            xml_entry.add_error_info(message=test_error.stderr)
            xml_parser.test_cases.append(xml_entry)
            continue

        other_out = ''
        for test_output in test_proc.stdout.splitlines():
            try:
                tests_results[test_output[:4]][test.name] += 1
            except KeyError:
                print(f'\t\t\t{test_output.lstrip()}')
                other_out += test_output.lstrip() + '\n'
                continue

            print(f'\t\t{test_output}')

            try:
                # Look for the test output message enclosed between brackets
                test_msg = re.search(r'(?<=\[).+(?=\])', test_output).group()
            except AttributeError:
                # If brackets are missing, extract output message
                # independently of an eventual error message
                test_msg = test_output[6:].split('; error:')[0]

            xml_entry = TestCase(test_msg, classname=test_id,
                                 elapsed_sec=wall_time,
                                 status=test_output[:4])

            try:
                fail_msg = re.search('(?<=error:).+', test_output).group()
                if not fail_msg.lstrip():
                    fail_msg += 'Empty message'
                if test_output[:4] == 'Warn':
                    xml_entry.add_failure_info(message=fail_msg.lstrip(),
                                               failure_type='warning')
                else:
                    xml_entry.add_failure_info(message=fail_msg.lstrip())
            except AttributeError:
                pass

            if test_proc.stderr:
                xml_entry.stderr = test_proc.stderr

            if other_out:
                xml_entry.stdout = other_out[:-1]
                other_out = ''

            xml_parser.test_cases.append(xml_entry)

    display_results(tests_results, error_list, skip_list)

    with open(xml_outfile, 'w') as fid:
        xml_parser.to_file(fid, [xml_parser])


def add_path_to_environment_variable(env_var, env_path):
    try:
        if str(env_path) not in os.environ[env_var]:
            os.environ[env_var] += f':{env_path}'
    except KeyError:
        os.environ[env_var] = str(env_path)


parser = argparse.ArgumentParser(description='Fluidity Unittest Harness')
parser.add_argument('-d', '--dir', default='bin/tests', metavar='',
                    help="""directory where to look for unittest binary(ies) -
                            defaults to '%(default)s'""")
parser.add_argument('-t', '--test', metavar='',
                    help="""single unittest to run - should match an existing
                            unittest within --dir""")
parser.add_argument('-f', '--from-file', metavar='',
                    help="""file containing unittests to run - expects one test
                            entry per line""")
parser.add_argument('-x', '--xml-output', metavar='',
                    help="""filename where to write the xml output - junit_xml
                            is required to generate the xml output""")
parser.add_argument('-c', '--clean', action='store_true',
                    help='remove the directory provided through --dir')
parser.add_argument('--efence', action='store_true',
                    help='links against libefence')
args = parser.parse_args()

fluidity_root = Path(sys.argv[0]).resolve().parent.parent
add_path_to_environment_variable('PYTHONPATH', fluidity_root / 'python')
add_path_to_environment_variable('PYTHONPATH', fluidity_root / 'lib')
if args.efence:
    add_path_to_environment_variable('LD_PRELOAD',
                                     Path('/usr/lib/libefence.so.0.0'))
# os.environ['EF_DISABLE_BANNER'] = '1'

if args.dir == 'bin/tests':
    tests_dir = fluidity_root / args.dir
else:
    tests_dir = Path(args.dir).resolve()
assert tests_dir.is_dir(), f'{tests_dir} is not a valid directory'

if args.clean:
    rmtree(tests_dir)
else:
    if args.test:
        tests = [tests_dir / args.test]
    elif args.from_file:
        with open(args.from_file, 'r') as fid:
            tests = [tests_dir / test_name.rstrip() for test_name in fid]
    else:
        tests = tests_dir.glob('test*')
    if args.xml_output:
        from junit_xml import TestSuite, TestCase
        if sys.version_info.major == 3 and sys.version_info.minor >= 8:
            from importlib.metadata import version
            assert float(version('junit_xml')) >= 1.9, (
                'Please update junit_xml')
        else:
            from pkg_resources import get_distribution
            assert float(get_distribution('junit_xml').version) >= 1.9, (
                'Please update junit_xml')
        unittest_harness(tests, args.xml_output)
    else:
        unittest_harness_no_output(tests)
