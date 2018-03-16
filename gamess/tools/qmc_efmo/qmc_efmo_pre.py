#!/usr/bin/env python

import sys, re

fname=sys.argv[1]
fname0=fname.split('.')[0]
fname=fname0+".dat"

key=re.compile(r'QMC CURRENT.*=\s*(\d+).*=\s*(\d+).*=\s*(\d+)')

f=open(fname)
lines=f.readlines()
f.close()

f=open("templ_files/head")
head_lines=f.readlines()
f.close()

f=open("templ_files/H.bas")
Hbas_lines=f.readlines()
f.close()

f=open("templ_files/O.bas")
Obas_lines=f.readlines()
f.close()

f=open("templ_files/H.ecp")
Hecp_lines=f.readlines()
f.close()

f=open("templ_files/O.ecp")
Oecp_lines=f.readlines()
f.close()

f_t1=open("temp1",'w')
f_t2=open("temp2",'w')

n=len(lines)
for ind in range(n):
   line=lines[ind]
   match=key.search(line)
   if match:
      msg=match.groups()
      i,j,k=map(int,[msg[0],msg[1],msg[2]])
      nat=int(lines[ind+3][:-1].split("=")[1])
      dname="wat_"+str(i)+"_"+str(j)+"_"+str(k)
      if j==0:
         print >>f_t1,dname
      else:
         print >>f_t2,dname
      fname=dname+".inp"
      f=open(fname,'w')
      for line in head_lines:
         print >>f,line[:-1]
      atoms=[]
      for ind2 in range(nat):
         pass
         v2=lines[ind+11+ind2].split()
         if j==0:
            a,c,x,y,z,e=v2[0][3:],float(v2[1]),float(v2[2]),float(v2[3]),float(v2[4]),float(v2[5])
            print >>f,a,c,x,y,z
            if a=="H":
               for line in Hbas_lines:
                  print >>f,line[:-1]
            if a=="O":
               for line in Obas_lines:
                  print >>f,line[:-1]
            atoms.append([a,c,x,y,z])
         else:
            a,c,x,y,z=v2[0][3:],float(v2[1]),float(v2[2]),float(v2[3]),float(v2[4])
            print >>f,a,c,x,y,z
            if a=="H":
               for line in Hbas_lines:
                  print >>f,line[:-1]
            if a=="O":
               for line in Obas_lines:
                  print >>f,line[:-1]
            atoms.append([a,c,x,y,z])
      print >>f," $END"
      print >>f," $ECP"
      H_flag=0
      O_flag=0
      for atom in atoms:
         [a,c,x,y,z]=atom
         if a=="H":
            if H_flag==0:
               for line in Hecp_lines:
                  print >>f,line[:-1]
               H_flag=1
            else:
               print >>f,"H-QMC GEN"
         elif a=="O":
            if O_flag==0:
               for line in Oecp_lines:
                  print >>f,line[:-1]
               O_flag=1
            else:
               print >>f,"O-QMC GEN"
      print >>f," $END"
      f.close()
f_t1.close()
f_t2.close()
