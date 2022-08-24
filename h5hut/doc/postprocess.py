#!/usr/bin/env python
import re

from lxml import etree
from lxml.html import builder as E

tree = etree.parse("tagfile")
tags = tree.xpath('//member[@kind="function"]')
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
    if re.search("^OUTPUT_DIRECTORY", line):
        prefix = line.split("=")[-1].strip()

if not prefix:
    prefix = "."

parser = etree.HTMLParser()
for filename, anchors in list(tags_per_files.items()):
    print("Processing file: {}/{}".format(prefix, filename))
    tree = etree.parse("{}/{}".format(prefix, filename), parser)
    for name, id in anchors:
        anchor = tree.xpath('//a[@id="%s"]' % id)
        el = E.A(id="%s" % name)
        anchor[0].addnext(el)

    f = open("{}/{}".format(prefix, filename), "w")
    f.write(etree.tostring(tree))
