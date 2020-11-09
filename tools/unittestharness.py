#!/usr/bin/env python3

import argparse
from collections import Counter
import os
from pathlib import Path
import re
from shutil import rmtree
import subprocess
import sys
import time

try:
    from junit_xml import TestSuite, TestCase
except ImportError:
    class TestSuite(object):
        def __init__(self, name, test_cases):
            self.test_cases = test_cases

        def to_file(self, *args):
            print('Cannot generate xml report without junit_xml module')

    class TestCase(object):
        def __init__(self, *args, **kwargs):
            pass

        def add_failure_info(self, *args, **kwargs):
            pass


def unitTestHarness(testsPath, unittest, from_file, xml_outfile):

    tests = []
    if unittest is not None:
        tests = [testsPath / unittest]
    elif from_file is not None:
        # Does not verify that the test exists. It only parses it to a list.
        # If the test does not exist, the unittestharness will simply skip it.
        try:
            with open(from_file, 'r') as f:
                tests = [testsPath / t for t in f.read().splitlines()]
        except IOError:
            sys.stderr.write(f'Could not read file: {from_file}\n')
            sys.exit(1)
    else:
        tests = testsPath.glob('test*')

    xmlParser = TestSuite('UnitTestHarness')
    # Pass, Warn and Fail are defined in femtools/Unittest_tools.F90
    # Make sure any change there is reflected here
    testsResults = {'Pass': Counter(), 'Warn': Counter(), 'Fail': Counter()}
    errorList = []

    for test in tests:
        print(f'\t{test.name}: Running')
        begin = time.process_time()
        testProc = subprocess.run(
            './' + test.name, shell=True, cwd=testsPath, encoding='utf-8',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        wallTime = time.process_time() - begin

        testId = f'unittest.{test.resolve().parts[-3]}.{test.name}'

        if testProc.stderr:
            if re.search(r'^\s*error[:\s]', testProc.stderr,
                         re.IGNORECASE | re.MULTILINE):
                print(f'ERROR: {test.name} failed')
                errorList.append(test.name)
                xmlEntry = TestCase('', classname=testId,
                                    elapsed_sec=wallTime, status='Error')
                xmlEntry.add_error_info(message=testProc.stderr)
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
                                elapsed_sec=wallTime, status=testOutput[:4])

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

    print(f'''RESULTS
\tPasses: {sum(testsResults['Pass'].values())}\n
\tWarnings: {sum(testsResults['Warn'].values())}
\t\t{list(testsResults['Warn'].keys())}\n
\tFailures: {sum(testsResults['Fail'].values())}
\t\t{list(testsResults['Fail'].keys())}\n
\tErrors: {len(errorList)}
\t\t{errorList}''')

    if xml_outfile is not None:
        with open(xml_outfile, 'w') as fid:
            xmlParser.to_file(fid, [xmlParser])
    return


parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''Fluidity UnitTestHarness''')
parser.add_argument('-d', '--dir', default='bin/tests', metavar='',
                    help="UnitTest binary location - defaults to 'bin/tests'")
parser.add_argument("-t", "--test", dest="unittest",
                    help="run a single unittest (by test name)", default=None)
parser.add_argument("--from-file", dest="from_file", default=None,
                    help="run tests listed in FROM_FILE (one test per line)")
parser.add_argument('-x', '--xml-output', metavar='',
                    help='Where to store the xml output')
parser.add_argument('-c', '--clean', action='store_true',
                    help='Remove UnitTests folder')
parser.add_argument('--efence', action='store_true',
                    help='Links against libefence')
args = parser.parse_args()

pyPath = (Path(sys.argv[0]).parent / '../python').resolve()
try:
    if str(pyPath) not in os.environ['PYTHONPATH']:
        os.environ['PYTHONPATH'] += f':{pyPath}'
except KeyError:
    os.environ['PYTHONPATH'] = str(pyPath)

ldPath = (pyPath / '../lib').resolve()
try:
    if str(ldPath) not in os.environ['LD_LIBRARY_PATH']:
        os.environ['LD_LIBRARY_PATH'] += f':{ldPath}'
except KeyError:
    os.environ['LD_LIBRARY_PATH'] = str(ldPath)

if args.efence:
    # Should be handled the same way as above
    os.environ['LD_PRELOAD'] = '/usr/lib/libefence.so.0.0'
    # os.environ['EF_DISABLE_BANNER'] = '1'

testsDir = pyPath / f'../{args.dir}'

if args.clean and testsDir.is_dir():
    rmtree(testsDir)
else:
    unitTestHarness(testsDir, args.unittest, args.from_file, args.xml_output)
