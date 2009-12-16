#!/usr/bin/env python

import sys
import os

svnurl = "http://amcg.ese.ic.ac.uk/svn/fluidity/trunk"
svnurls = {}
svnurls['libadapt']      = "http://amcg.ese.ic.ac.uk/svn/libadapt/trunk"
svnurls['libadaptivity'] = "http://amcg.ese.ic.ac.uk/svn/libadaptivity/trunk"
svnurls['libsam']        = "http://amcg.ese.ic.ac.uk/svn/libsam/trunk"
svnurls['libvtkfortran'] = "http://amcg.ese.ic.ac.uk/svn/vtkfortran/trunk"

def get_svn_head():
  f = os.popen("svn info " + svnurl + " | grep Revision | awk '{print $NF}'")
  return int(f.read())

def get_svn_author(revision):
  f = os.popen("svn info -r " + str(revision) + " " + svnurl + "| grep 'Last Changed Author' | awk '{print $NF}'")
  return f.read()[:-1]

def get_svn_date(revision):
  f = os.popen("svn info -r " + str(revision) + " " + svnurl + "| grep 'Last Changed Date'")
  line = f.read()
  return '{' + ' '.join(line.split()[3:6]) + '}'

def set_up_repo(file, dir, revision):
  print "Setting up repository: revision %s" % p
  os.system("rm -rf fluidity-" + str(revision))
  f = os.popen("svn co -r " + str(revision) + " " + svnurl + " fluidity-" + str(revision))
  for line in f:
    print line,
  os.system("mkdir -p fluidity-" + str(revision) + os.sep + "tests")
  os.system("cp " + dir + "tools/*.py fluidity-" + str(revision) + os.sep + "tools") # copy the test harness
  os.system("rm -rf fluidity-" + str(revision) + os.sep + "tests" + os.sep + file[:-4])
  os.system("cp -r " + dir + "tests" + os.sep + file[:-4] + " fluidity-" + str(revision) + os.sep + "tests")
  os.chdir("fluidity-" + str(revision))
  date = get_svn_date(revision)
  for lib in ['libadapt', 'libadaptivity', 'libsam', 'libvtkfortran']:
    if(os.path.exists(lib)):
      os.chdir(lib)
      os.system("svn merge -r 'HEAD:%s' %s" % (date, svnurls[lib]))
      os.chdir(os.pardir)
  os.system("sed -i 's@./lib/libadapt @./lib/libadapt.a @' configure.in ; autoconf")
  os.chdir(os.pardir)
  print

if __name__ == "__main__":
  import optparse

  parser = optparse.OptionParser()
  parser.add_option("-f", "--file", dest="file", help="specific test case to run (by filename)", default="")
  parser.add_option("-r", "--revision", dest="revision", help="last known good revision", default="")
  parser.add_option("-b", "--broken", dest="broken", help="known broken revision", default=get_svn_head())

  (options, args) = parser.parse_args()

  if not hasattr(options, "file") or not hasattr(options, "revision"):
    print "Error: look at the help."
    sys.exit(1)

  low = int(options.revision)
  high = int(options.broken)
  file = options.file
  dir = os.getcwd().split(os.sep)[-2] + os.sep
  os.chdir(os.pardir + os.sep + os.pardir)
  low_last_changed = False

  print "Performing a binary search on %s: (low, high) == (%s,%s)" % (file, low, high)

  p = (low + high) / 2
  while (low <= high):
    set_up_repo(file, dir, p)
    os.chdir("fluidity-" + str(p))
    maxret = 0
    ret = os.system("./configure --enable-2d-adaptivity"); maxret = max(ret, 0)
    ret = os.system("make"); maxret = max(ret, 0)
    os.chdir("tools")
    os.chdir(os.pardir)
    if (maxret != 0):
      print "WARNING: revision %s failed to compile" % p
      os.chdir(os.pardir)
      if low_last_changed: p = p - 1
      else: p = p + 1
      continue
    ret = os.system("tools/testharness.py -f " + file)
    os.chdir(os.pardir)
    if (ret == 0):
      low = p + 1
      low_last_changed = True
    else: 
      high = p - 1
      low_last_changed = False
    p = (low + high) / 2

  if ret == 0:
    p = p + 1

  print "The test case broke on revision: %s" % (p)
  print "Author: %s" % get_svn_author(p) 
