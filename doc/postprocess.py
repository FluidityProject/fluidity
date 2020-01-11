#!/usr/bin/env python

from lxml import etree
from lxml.html import builder as E
import os
import re

tree   = etree.parse('tagfile')
tags = tree.xpath ('//member[@kind="function"]')
tags_per_files = {}
for tag in tags:
     name = tag.find("name").text
     anchor = tag.find("anchor").text
     anchorfile = tag.find("anchorfile").text
     if anchorfile in tags_per_files:
         tags_per_files[anchorfile].append([name, anchor])
     else:
         tags_per_files[anchorfile] = [[name, anchor]]


prefix = None
for line in open("Doxyfile", "r"):
     if re.search ("^OUTPUT_DIRECTORY", line):
          prefix = line.split("=")[-1].strip()

if not prefix:
     prefix = "."

parser = etree.HTMLParser()
for filename,anchors in tags_per_files.items():
     print "Processing file: %s/%s" % (prefix, filename)
     tree = etree.parse ("%s/%s" % (prefix, filename), parser)
     for name,id in anchors:
          anchor = tree.xpath ('//a[@id="%s"]' % id)
          el = E.A(id="%s" % name)
          anchor[0].addnext (el)

     f = open ("%s/%s" % (prefix, filename), "w")
     f.write (etree.tostring (tree))
