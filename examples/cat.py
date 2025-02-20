#!/usr/bin/python
# cat with blank lines
# mcopy -o -v cat.py a:

from sys import argv

for fn in argv[1:]:
  f=open(fn,'r')
  for line in f:
    print line[:-1]
  print
  f.close()
