import os
import sys
import yaml

from   findfile import file_check, findfile, fetch_header, supported
from   config   import smart_open, CustomDumper

def diagnose():
   result = fetch_header(allsupported=True)

   for key in result.keys():
      with smart_open() as fh:
         yaml.dump(result[key], fh, default_flow_style=False, sort_keys=False)
         fh.write('\n')

def summary(fpath=None):
    if not fpath:
        fpath = findfile('summary_log')
    
    smart_open(fpath)
    
    diagnose()

    file_check()  

def tidyup():
    summary()

if __name__ == '__main__':
   tidyup()

   print('\n\nDone.\n\n')
