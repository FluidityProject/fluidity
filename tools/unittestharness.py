#!/usr/bin/env python3

# Pass, Warn and Fail are defined in femtools/Unittest_tools.F90
# Make sure any changes there are reflected here

import argparse
from collections import Counter
import os
from pathlib import Path
import re
from shutil import rmtree
import subprocess
import sys
from time import process_time


def printResults(testsResults, errorList, skipList):
    print(f'''
###########
# Results #
###########

\tPasses: {sum(testsResults['Pass'].values())}\n
\tWarnings: {sum(testsResults['Warn'].values())}
\t\t{list(testsResults['Warn'].keys())}\n
\tFailures: {sum(testsResults['Fail'].values())}
\t\t{list(testsResults['Fail'].keys())}\n
\tErrors: {len(errorList)}
\t\t{errorList}\n
\tSkipped: {len(errorList)}
\t\t{skipList}''')


def unitTestHarnessNO(tests):
    testsResults = {'Pass': Counter(), 'Warn': Counter(), 'Fail': Counter()}
    errorList, skipList = [], []

    for test in tests:
        print(f'\t{test.name}: Running')
        testProc = subprocess.run(
            './' + test.name, shell=True, cwd=test.parent, encoding='utf-8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if re.search(r'^\s*error[:\s]', testProc.stderr,
                     re.IGNORECASE | re.MULTILINE):
            print(f'ERROR: {test.name} failed')
            errorList.append(test.name)
            continue
        elif ': not found' in testProc.stderr:
            print(f'WARNING: {test.name} not found')
            skipList.append(test.name)
            continue

        for testOutput in testProc.stdout.splitlines():
            try:
                testsResults[testOutput[:4]][test.name] += 1
            except KeyError:
                print(testOutput.lstrip())
                continue

            print(f'\t\t{testOutput}')

    printResults(testsResults, errorList, skipList)


def unitTestHarness(tests, xml_outfile):
    xmlParser = TestSuite('UnitTestHarness')
    testsResults = {'Pass': Counter(), 'Warn': Counter(), 'Fail': Counter()}
    errorList, skipList = [], []

    for test in tests:
        print(f'\t{test.name}: Running')
        begin = process_time()
        testProc = subprocess.run(
            './' + test.name, shell=True, cwd=test.parent, encoding='utf-8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        wallTime = process_time() - begin

        testId = f'unittest.{test.resolve().parts[-3]}.{test.name}'

        if re.search(r'^\s*error[:\s]', testProc.stderr,
                     re.IGNORECASE | re.MULTILINE):
            print(f'ERROR: {test.name} failed')
            errorList.append(test.name)
            xmlEntry = TestCase('', classname=testId,
                                elapsed_sec=wallTime, status='Error')
            xmlEntry.add_error_info(message=testProc.stderr)
            xmlParser.test_cases.append(xmlEntry)
            continue
        elif ': not found' in testProc.stderr:
            print(f'WARNING: {test.name} not found')
            skipList.append(test.name)
            xmlEntry = TestCase('', classname=test.name,
                                elapsed_sec=wallTime, status='Skipped')
            xmlEntry.add_skipped_info(message='Not found')
            xmlParser.test_cases.append(xmlEntry)
            continue

        otherOut = ''
        for testOutput in testProc.stdout.splitlines():
            try:
                testsResults[testOutput[:4]][test.name] += 1
            except KeyError:
                print(testOutput.lstrip())
                otherOut += testOutput.lstrip() + '\n'
                continue

            print(f'\t\t{testOutput}')

            try:
                testMsg = re.search(r'(?<=\[).+(?=\])', testOutput).group()
            except AttributeError:
                testMsg = testOutput[6:].split('; error:')[0]

            xmlEntry = TestCase(testMsg, classname=testId,
                                elapsed_sec=wallTime,
                                status=testOutput[:4])

            try:
                failMsg = re.search('(?<=error:).+', testOutput).group()
                if not failMsg.lstrip():
                    failMsg += 'Empty message'
                if testOutput[:4] == 'Warn':
                    xmlEntry.add_failure_info(message=failMsg.lstrip(),
                                              failure_type='warning')
                else:
                    xmlEntry.add_failure_info(message=failMsg.lstrip())
            except AttributeError:
                pass

            if testProc.stderr:
                xmlEntry.add_error_info(message=testProc.stderr,
                                        error_type='warning')

            if otherOut:
                xmlEntry.stdout = otherOut[:-1]
                otherOut = ''

            xmlParser.test_cases.append(xmlEntry)

    printResults(testsResults, errorList, skipList)

    with open(xml_outfile, 'w') as fid:
        xmlParser.to_file(fid, [xmlParser])


def updateEnvVar(envVar, envPath):
    try:
        if str(envPath) not in os.environ[envVar]:
            os.environ[envVar] += f':{envPath}'
    except KeyError:
        os.environ[envVar] = str(envPath)


parser = argparse.ArgumentParser(description='''Fluidity UnitTestHarness''')
parser.add_argument('-d', '--dir', default='bin/tests', metavar='',
                    help="""directory where to look for UnitTest binary(ies) -
                            defaults to '%(default)s'""")
parser.add_argument('-t', '--test', metavar='',
                    help='''single UnitTest to run - should match an existing
                            UnitTest within --dir''')
parser.add_argument('-f', '--from-file', metavar='',
                    help='''file containing UnitTests to run - expects one test
                            entry per line''')
parser.add_argument('-x', '--xml-output', metavar='',
                    help='''filename where to write the xml output - junit_xml
                            is required to generate the xml output''')
parser.add_argument('-c', '--clean', action='store_true',
                    help='remove the directory provided through --dir')
parser.add_argument('--efence', action='store_true',
                    help='links against libefence')
args = parser.parse_args()

pyPath = (Path(sys.argv[0]).parent / '../python').resolve()
updateEnvVar('PYTHONPATH', pyPath)
updateEnvVar('LD_LIBRARY_PATH', (pyPath / '../lib').resolve())
if args.efence:
    updateEnvVar('LD_PRELOAD', Path('/usr/lib/libefence.so.0.0'))
# os.environ['EF_DISABLE_BANNER'] = '1'

testsDir = pyPath / f'../{args.dir}'
assert testsDir.is_dir(), f'{testsDir.resolve()} is not a valid directory'

if args.clean:
    rmtree(testsDir)
else:
    if args.test:
        tests = [testsDir / args.test]
    elif args.from_file:
        with open(args.from_file, 'r') as fid:
            tests = [testsDir / testName.rstrip() for testName in fid]
    else:
        tests = testsDir.glob('test*')
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
        unitTestHarness(tests, args.xml_output)
    else:
        unitTestHarnessNO(tests)
