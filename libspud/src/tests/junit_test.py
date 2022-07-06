#!/usr/bin/env python3
import glob
import os
import subprocess

from junit_xml import TestCase
from junit_xml import TestSuite

tests = glob.glob("bin/*")

suites = []

for test in tests:
    output = subprocess.check_output(os.path.abspath(test)).decode()
    title = ""
    for line in output.split("\n"):
        if not line.strip():
            continue
        if line[1:4] == "***":
            title = line.replace("***", "").strip()
            suites.append(TestSuite("%s" % title, []))
        else:
            result, name = [part.strip() for part in line.split(":", 1)]
            result = result == "Pass"
            if not result:
                name, error = name.split(";")
            name = name[1:-1]
            suites[-1].test_cases.append(TestCase(name, "{}.{}".format(title, name)))
            if not result:
                suites[-1].test_cases[-1].add_failure_info("Fail", error)

with open("test_results.xml", "w") as handle:
    suites[-1].to_file(handle, suites)
