#!/usr/bin/env python

import re
from subprocess import PIPE, Popen

def cmdline(command):
    process = Popen(
        args=command,
        stdout=PIPE,
        shell=True
    )
    return process.communicate()[0]

def emin(jname):
   cmd="/global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/nexus/executables/qmca -q e "+jname+".s*.scalar.dat"
   out=cmdline(cmd)
   olines=out[2:-1].split("\n")
   n_min=0
   e_min=0.
   for oline in olines:
      n,e=int(oline.split()[2]),float(oline.split()[5])
      #print n,e
      if e<e_min: 
         n_min=n 
         e_min=e
   n_min=n_min-1
   if n_min<10:
      return "s00"+str(n_min)
   else:
      return "s0"+str(n_min)

def bscript(fname,jname,nproc,time,cmd):
   f=open(fname,'w')
   print >>f,"#!/bin/bash -l"
   print >>f,"#SBATCH -J "+jname
   print >>f,"#SBATCH -p regular"
   print >>f,"#SBATCH -N "+str(nproc)
   print >>f,"#SBATCH -t "+time
   print >>f,"#SBATCH -A m2862"
   print >>f,"#SBATCH --qos=premium"
   print >>f
   print >>f,"cd $SLURM_SUBMIT_DIR"
   print >>f
   print >>f, cmd
   f.close()

def oinp(fname,jname,wname,templl):
   f=open(fname,'w')
   templl[2]='  <project id="'+jname+'" series="0"/>\n'
   templl[6]='  <include href="'+wname+'"/>\n'
   for line in templl:
      print >>f,line[:-1]

def addJ3(fname,J3_templl):
   f=open(fname,'r')
   lines=f.readlines()
   f.close()
   f=open(fname,'w')
   for line in lines[:-2]:
      print >>f,line[:-1]
   for line in J3_templl:
      print >>f,line[:-1]
   for line in lines[-2:]:
      print >>f,line[:-1]
   f.close()

def dmc_e(dname):
   #cmd="/global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/nexus/executables/qmca -q e "+dname+"/ts-h2o.s001.scalar.dat"
   cmd="/global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/nexus/executables/qmca -q e ts-h2o.s001.scalar.dat"
   out=cmdline(cmd)
   olines=out.split()
   return float(olines[5]),float(olines[7])

def gamess_e(fname):
   key=re.compile(r'TOTAL ENERGY =\s*(\S+)')
   f=open(fname,'r')
   lines=f.readlines()
   f.close()
   for line in lines:
      match=key.search(line)
      if match:
         return float(match.groups()[0])
