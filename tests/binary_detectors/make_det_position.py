#!/usr/bin/env python3

import struct

""" 
Create initial positions for detectors in binary to test the from_file functionality
"""

f = open('positions.dat', 'wb')
f.write(struct.pack('d',0.0))
f.write(struct.pack('d',0.0))
f.write(struct.pack('d',0.25))
f.write(struct.pack('d',0.25))
f.write(struct.pack('d',0.4))
f.write(struct.pack('d',0.4))
f.close()
