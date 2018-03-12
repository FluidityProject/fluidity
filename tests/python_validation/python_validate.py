#/usr/bin/env python2

""" 

python2 python_validate.py [-p] [-f] [file1 file2 file3]

or 

python3 python_validate.py [-p] [-f] [file1 file2 file3]

-f  (--fixup)  Automatically attempt to fix obvious python3 errors using reindent and 2to3.
-p  (--pep8)   Run flake8 to test adherence to pep8 coding stanaards.

Script to attempt to validate python files and code snippets in the Fluidity code base against python2 and python3 syntax.

"""

from __future__ import print_function
import sys
import os
import ast
import subprocess
import tempfile

from lxml import etree

import argparse

def flake8_snippet(name, snippet):
    """Run flake8 against the snippet passed in, with reporting name as 'name'."""
    with tempfile.TemporaryFile(mode='w') as temp: 
        temp.write(snippet)
        temp.seek(0)
        try:
            subprocess.check_output(['flake8 --stdin-display-name="%s" -'%name],
                                    stdin = temp, shell=True)
        except subprocess.CalledProcessError as err:
            print(err.output.decode("utf8"))

def print_as_function(filename, add_future_import=False):
    """Handrolled function to turn
    print blah
       into
    print(blah)"""
    lines = []

    def sanitize(line, sep):
        if len(sep)==0:
            if 'print ' in line:
                return line.replace('print ', 'print(')+')'
            else:
                return line
        else:
            return sep[0].join([sanitize(_, sep[1:]) for _ in line.split(sep[0])])

    with open(filename, 'r') as temp:
        for line in temp.readlines():
            if len(line.strip())==0 or line.strip()[0]=='#':
                pass
            if 'print ' in line:
                line = sanitize(line.rstrip('\n'), [';', '#'])+'\n'
            lines.append(line)
    with open(filename, 'w') as temp:
        if add_future_import:
            temp.write('from __future__ import print_function')
        temp.writelines(lines)

def fix_snippet(name, snippet):
    temp, temp_name = tempfile.mkstemp()
    os.write(temp, snippet.encode("utf8"))
    os.close(temp)
    #use the local reindent script to update to 4 spaces and no tabs
    subprocess.call(['python2 -m reindent -n %s'%temp_name],
                    shell=True)
    #run 2to3 to 
    print_as_function(temp_name)
    subprocess.call(['2to3 -x future -npw %s'%temp_name],
                        shell=True)
    with open(temp_name, 'r') as temp:
        out = temp.read()
    os.remove(temp_name)
    if out == snippet:
        out == None
    return out
    

def get_name(element):
    return "/"+"/".join(reversed([element.tag]+[s.tag for s in element.iterancestors()]))

def process_snippet(name, snippet):
    """Parse and validate a string containing a snippet of Python code.

    Return a Boolean.

    On failure, also print the string with line numbers.
    """ 
    fail=False
    snippet = snippet.strip()
    try:
        ast.parse(snippet)
    except SyntaxError as err:
        print(str(err))
        fail = True
        lines = snippet.split('\n')
        nchar = len(str(len(lines)))
        for k, line in enumerate(lines):
            print(("%{0}d:".format(nchar))%(k+1), repr(line)[1:-1])
    return fail

def validate_file(filename, fix):
    """Run the validotor on a whole file."""
    with open(filename, 'r') as testfile:
        failed = process_snippet(filename, testfile.read())
        if failed:
            if fix:
                print_as_function(filename)
            return True
    return False


def validate_xml_snippets(filename, test_pep8, fix):
    """Run the validotor on elements with attribute 'language=python' in an xml file."""

    snippets=[]

    with open(filename, 'r') as xmlfile:
        tree = etree.ElementTree(file=filename)


        for ele in tree.iter():
            if ele.attrib.get('language',None)=='python':
                snippets.append((get_name(ele), ele.text, ele))

        failed = False
        modified = False
        for name, snippet, ele in snippets:
            if(test_pep8):
                flake8_snippet(name, snippet)
            if process_snippet(name, snippet):
                failed = True
                if fix:
                    fixed = fix_snippet(name, snippet)
                    if fixed:
                        modified = True
                        ele.text = fixed
                else:
                    return True
        if failed and modified:
            tree.write(filename, xml_declaration=True, encoding='utf-8',
                   with_tail=False, pretty_print=True)

    return failed

def main():
    """Main python validation routine."""

    parser = argparse.ArgumentParser(usage="<python_executable> python_validate.py [-p] [-f] [file1 file2 file3]")
    parser.add_argument("-f", "--fixup",  action="store_true", dest="fix", default=False,
                  help="Automatically attempt to fix obvious python3 errors using reindent and 2to3")
    parser.add_argument("-p", "--pep8", action="store_true", dest="test_pep8", default=False,
                  help="Run flake8 to test adherence to pep8 coding stanaards.")
    parser.add_argument('files', nargs='*',
                    help='files to process')
    options = parser.parse_args()
    failures = []

    failed = False
    print(options.files)
    for filename in options.files:
        if filename.rsplit('.',1)[-1]=='py':
            if validate_file(filename, options.fix):
                failed =True
                failures.append(filename)
        else:
            if validate_xml_snippets(filename, 
                                     test_pep8=options.test_pep8,
                                     fix=options.fix):
                failed =True
                failures.append(filename)

    if failed:
        print("Failed:")
        for filename in failures:
            print(filename)
        exit(1)
    else:
        print("Passed!")
        exit(0)

def check_versions():
    versions = []
    for ver in ("python2", "python3"):
        try:
            subprocess.call([ver+' --version'], shell=True)
            versions.append(ver)
        except subprocess.CalledProcessError:
            pass
    return versions
    

if __name__ == "__main__":
    main()



