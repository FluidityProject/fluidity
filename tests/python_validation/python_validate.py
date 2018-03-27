#/usr/bin/env python2

""" 

python2 python_validate.py [-p] [-f] [file1 file2 file3]

or 

python3 python_validate.py [-p] [-f] [file1 file2 file3]

-f  (--fixup)  Automatically attempt to fix obvious python3 errors using reindent and 2to3.
-p  (--pep8)   Run flake8 to test adherence to pep8 coding standards.

Script to attempt to validate python files and code snippets in the Fluidity code base against python2 and python3 syntax.

"""

from __future__ import print_function
import sys
import os
import ast
import subprocess
import pickle
import io
import tokenize
import tempfile

from xml.etree.ElementTree import ElementTree

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

    def isplit(line, indices):
        return [line[i+1:j] for i, j in zip(indices, indices[1:]+[None])]

    def split_comments(line):
        tokens = tokenize.generate_tokens(io.StringIO(line).readline)
        indices = [-1]
        for _ in tokens:
            if _[0] == tokenize.COMMENT:
                return isplit(line.replace(_[1],''), indices), _[1]
            if _[0] == tokenize.OP and _[1] == ';':
                indices.append(_[2][1])
                continue
        return isplit(line, indices), ''

    def sanitize(line):
        if 'print ' in line:
            return line.replace('print ', 'print(')+')'
        else:
            return line

    with open(filename, 'r') as temp:
        for line in temp.readlines():
            if len(line.strip())==0 or line.strip()[0]=='#':
                pass
            if 'print ' in line:
                line, comment = split_comments(line)
                line = ';'.join([sanitize(_.strip('\n')) for _ in line])+comment+'\n'
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
    fail = False
    snippet = snippet.strip()
    try:
        ast.parse(snippet)
    except SyntaxError as err:
        fail = "%s\n%s\n"%(name, err)
        lines = snippet.split('\n')
        nchar = len(str(len(lines)))
        for k, line in enumerate(lines):
            fail += ("%{0}d: %s\n".format(nchar))%((k+1), repr(line)[1:-1])
        print(fail)
    return fail

def validate_file(filename, fix):
    """Run the validator on a whole file."""
    with open(filename, 'r') as testfile:
        failed = process_snippet(filename, testfile.read())
        if failed:
            if fix:
                print_as_function(filename)
    return failed

def validate_xml_snippets(filename, test_pep8, fix):
    """Run the validator on elements with attribute 'language=python' in an xml file."""

    snippets=[]

    with open(filename, 'r') as xmlfile:
        tree = ElementTree(file=filename)

        #build iterator by hand, since not using lxml
        exp_tree = {}

        def expand_tree(ele, path):
            epath = path+ele.tag+'/'
            exp_tree[epath] = ele
            for _ in ele.getchildren():
                expand_tree(_, epath)

        expand_tree(tree.getroot(), '/')

        for name, ele in exp_tree.items():
            if ele.attrib.get('language',None)=='python':
                snippets.append((name, ele.text, ele))

        failed = ''
        modified = False
        for name, snippet, ele in snippets:
            name = name+ele.attrib.get('name','')
            if(test_pep8):
                flake8_snippet(name, snippet)
            test = process_snippet(name, snippet)
            if test:
                failed += test
                if fix:
                    fixed = fix_snippet(name, snippet)
                    if fixed:
                        modified = True
                        ele.text = fixed
                else:
                    return failed
        if failed and modified:
            tree.write(filename, xml_declaration=True, encoding='utf-8')

    return failed or False

def main():
    """Main python validation routine."""

    parser = argparse.ArgumentParser(usage="<python_executable> python_validate.py [-f] [-p] [-w output] [file1 file2 file3]")
    parser.add_argument("-f", "--fixup",  action="store_true", dest="fix", default=False,
                  help="Automatically attempt to fix obvious python3 errors using reindent and 2to3")
    parser.add_argument("-w", "--write-failures",  action="store", dest="write", default=None,
                  help="Write failures to disk.")
    parser.add_argument("-p", "--pep8", action="store_true", dest="test_pep8", default=False,
                  help="Run flake8 to test adherence to pep8 coding stanaards.")
    parser.add_argument('files', nargs='*',
                    help='files to process')
    options = parser.parse_args()
    failures = {}

#    print(options.files)
    for filename in options.files:
        print(filename)
        if filename.rsplit('.',1)[-1]=='py':
            test = validate_file(filename, options.fix)
            if test:
                failures[filename] = test
        else:
            test =  validate_xml_snippets(filename, 
                                     test_pep8=options.test_pep8,
                                     fix=options.fix)
            if test:
                failures[filename] = test

    if options.write:
        with open(options.write.strip(), 'wb') as output:
            pickle.dump(failures, output, protocol=2)

    if failures:
        print("Failed:")
        for filename in failures:
            print(filename)
        exit(100)
    else:
        print("Passed!")
        exit(0)

def check_versions():
    versions = []
    for ver in ("python2", "python3"):
        try:
            subprocess.check_output([ver+' --version'], shell=True)
            versions.append(ver)
        except subprocess.CalledProcessError:
            pass
    return versions
    

if __name__ == "__main__":
    main()



