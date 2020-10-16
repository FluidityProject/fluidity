#!/usr/bin/env python3

import os
import re
import glob
import time
import shutil

try:
  from junit_xml import TestSuite, TestCase
except ImportError:
    class TestSuite(object):
        def __init__(self, name, test_cases):
            self.test_cases=test_cases
        def to_file(self,*args):
            print("cannot generate xml report without junit_xml module.")
    class TestCase(object):
        def __init__(self,*args,**kwargs):
            pass
        def add_failure_info(self,*args,**kwargs):
            pass

class UnitTest:
    def __init__(self, exe):
        self.exe = os.curdir + os.sep + exe.split(os.sep)[-1]
        self.verbose = True
        self.dir = os.sep.join((os.curdir + os.sep + exe).split(os.sep)[:-1])
        self.cwd = os.getcwd()
        self.fail_msgs = []
        self.pass_msgs = []
        self.warn_msgs = []

    def log(self, msg):
        if self.verbose and msg != '':
            print("    %s: %s" % (self.exe.split('/')[-1], msg))

    def run(self):
        os.chdir(self.dir)
        # self.log("chdir " + self.dir)
        self.log("Running")
        f = os.popen(self.exe)
        self.output = f.read()
        os.chdir(self.cwd)
        
        exitStatus = f.close()
        if exitStatus is None:
          return 0
        else:
          return exitStatus

    def parse(self):
        passcount = 0
        warncount = 0
        failcount = 0

        for line in self.output.split('\n'):
            line = line.lstrip()
            self.log(line)

            # The junit integration is closely dependent to the output from
            # report_test in the Fortran Unitest_tools.F90 module, so if logging
            # in the Fortran module is ever altered remember to update these
            # splits as well accordingly.
            # Get the output after Pass, Warn, Fail
            if line.startswith("Pass"):
                self.pass_msgs.append(line.split("Pass: ", 1)[1])
                passcount = passcount + 1

            elif line.startswith("Warn"):
                self.warn_msgs.append(line.split("Warn: ", 1)[1])
                warncount = warncount + 1

            elif line.startswith("Fail"):
                self.fail_msgs.append(line.split("Fail: ", 1)[1])
                failcount = failcount + 1

        return (passcount, warncount, failcount)

class UnitTestHarness:
    def __init__(self, dir, xml_outfile=''):
        self.tests = []
        self.xml_parser = TestSuite('TestHarness',[])
        self.xml_outfile = xml_outfile
        self.xml_reports = []
        self.cwd = os.getcwd()
        self.dir = dir

        if dir[-1] == '/': dir = dir + "*"
        else: dir = dir + "/*"

        # Ignore all .vtu files in the test directory
        files = glob.glob(f'{dir}[!.vtu]')
        for file in files:
            if not os.path.isdir(file):
                self.tests.append(UnitTest(file))

    def clean(self):
        shutil.rmtree(os.path.join(self.cwd, 'bin/tests'))

    def run(self):
        passcount = 0
        warncount = 0
        failcount = 0

        warntests = []
        failtests = []
        passtests = []

        for test in self.tests:
            wall_time = time.time()
            exitStatus = test.run()
            
            (P, W, F) = test.parse()

            # Discard relative path characters
            test_name = test.exe[2:]

            # Get the location of the test
            # If it is a symlink get its original location in a relative path
            root_test_dir = ''
            t_path = os.path.join(self.cwd, self.dir, test_name)
            if os.path.islink(t_path):
                t_path = os.readlink(t_path)

            # This relies on the unittests being placed in a tests/ directory
            root_test_dir = ''
            # If symlink, path is relative, clean and extract substring between
            # the last ./ and /tests/
            try:
                root_test_dir = re.search('.+/(.+?)/tests/', t_path).group(1)
            # Not a symlink, so try and derive from binary dir location
            except AttributeError:
                root_test_dir = re.search('(.+?)/tests/', self.dir).group(1)
            except:
                root_test_dir = 'unknown-test-root'
            classname = f'unittest.{root_test_dir}.{test_name}'

            if (P, W, F) == (0, 0, 0):
                msg = "    WARNING: no output from test"
                print(msg)
                test.warn_msgs.append(msg)
                warncount += 1
                warntests.append(test_name)

            if P > 0:
                passtests.append(test_name)

                for t in test.pass_msgs:
                    tc = TestCase(t, classname, time.time() - wall_time)
                    self.xml_reports.append(tc)
                    self.xml_parser.test_cases += [tc]

            if W > 0:
                warntests.append(test_name)

                for t in test.warn_msgs:
                    # Extract the test name from the error message
                    name = t.split('; ', 1)[0]
                    msg = t.split('error: ', 1)[1]
                    tc = TestCase(name, classname, time.time() - wall_time)
                    tc.add_failure_info('Failure', msg)
                    self.xml_reports.append(tc)
                    self.xml_parser.test_cases += [tc]

            if F > 0:
                failtests.append(test_name)

                for t in test.fail_msgs:
                    # Extract the test name from the error message
                    name = t.split('; ', 1)[0]
                    msg = t.split('error: ', 1)[1]
                    tc = TestCase(name, classname, time.time() - wall_time)
                    tc.add_failure_info('Failure', msg)
                    self.xml_reports.append(tc)
                    self.xml_parser.test_cases += [tc]
              
            if not exitStatus == 0:
                print('    ERROR: non-zero exit code from test')
                failcount += 1
                if not test_name in failtests:
                    failtests.append(test_name)

            passcount += P
            warncount += W
            failcount += F

        print("RESULTS")
        print("    Passes:   %d" % passcount)
        if len(warntests) == 0:
            print("    Warnings: %d" % warncount)
        else:
            print("    Warnings: %d; tests = %s" % (warncount, warntests))
        if len(failtests) == 0:
            print("    Failures: %d" % failcount)
        else:
            print("    Failures: %d; tests = %s" % (failcount, failtests))

        if self.xml_outfile != '':
            fd = open(os.path.join(self.cwd, self.xml_outfile), 'w')
            self.xml_parser.to_file(fd, [self.xml_parser])
            fd.close()

if __name__ == "__main__":
    import sys
    import optparse

    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir', default='bin/tests', help='unittest binary location default bin/tests')
    parser.add_option('-x', '--xml-output', dest='xml_outfile', default='', help='filename for xml output')
    parser.add_option('-c', '--clean', action='store_true', dest='clean', default = False)
    parser.add_option('--electricfence', action='store_true', dest='electricfence', default = False, help='link against electric fence lib')
    (options, args) = parser.parse_args()

    try:
        os.environ["PYTHONPATH"] = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "..", "python")) + ":" + os.environ["PYTHONPATH"]
    except KeyError:
        os.putenv("PYTHONPATH", os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "..", "python")))

    try:
        os.environ["LD_LIBRARY_PATH"] = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "..", "lib")) + ":" + os.environ["LD_LIBRARY_PATH"]
    except KeyError:
        os.putenv("LD_LIBRARY_PATH", os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "..", "lib")))

    if options.electricfence:
        os.putenv("LD_PRELOAD", "/usr/lib/libefence.so.0.0")
        #os.putenv("EF_DISABLE_BANNER", "1")

    unittestharness = UnitTestHarness(dir=options.dir, 
                                      xml_outfile=options.xml_outfile)

    if options.clean:
        unittestharness.clean()
    else:
        unittestharness.run()
