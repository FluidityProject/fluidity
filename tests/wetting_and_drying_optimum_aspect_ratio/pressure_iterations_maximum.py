#!/usr/bin/env python

import numpy as np

def find_iterations_pressure(filename):
  f = open(filename, 'r')
  pit = []
  for line in f:
    if r'DeltaP PETSc n/o iterations:' in line:
      columns = line.strip().split()
      pit.append(int(columns[-1]))

  if len(pit) == 0:
    return None

  data = np.array(pit)
  entries = data.shape[0]
  totals = np.sum(data, axis=0) 
  maximums = np.max(data, axis=0)
  averages = totals / entries


  result = {
    'entries': entries,
    'total': totals,
    'maximum': maximums,
    'average': averages,
  }
  return result

def main():
  import sys
  args = sys.argv[1:]

  if len(args) == 0:
    print 'Please provide a filename.'
    sys.exit(1)

  filename = sys.argv[1]
  result = find_iterations_pressure(filename)


  if result is None:
    print 'No iteration counts found'
    sys.exit(1)

  print result['entries'], result['total'], result['maximum'], result['average']


if __name__ == "__main__":
  main()


#   Run,       a    ts_end      t_end    iua    ipa    ium    ipm
#     1,    0.01         1         0s    266    300    277    304
#     2,     1.0         1         0s    267     81    277     85
#     3,   100.0         1         0s    285    275    294    280

