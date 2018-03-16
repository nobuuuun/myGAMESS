#!/usr/bin/env python

from basic import dmc_e,gamess_e
import sys

dname=sys.argv[1]
   
en,er=dmc_e(dname) 
en0=gamess_e(dname)
print en0,en,er
