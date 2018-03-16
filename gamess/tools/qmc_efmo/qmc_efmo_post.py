#!/usr/bin/env python

import numpy as np

#e0=-69.068227
#err0=0.000053

#e0V=-69.007713
#err0V=0.000625 

H2kcal_mol=627.509

#b1
e0H=-1216.3414728073
e0eH=-1216.3420251981
#b2
#e0H=-1216.9913347491
#e0eH=-1216.9961177501

e1=np.zeros((16,))
e1_0=np.zeros((16,))
e1c=np.zeros((16,))
err1=np.zeros((16,))
e2=np.zeros((16,16))
e2_0=np.zeros((16,16))
e2c=np.zeros((16,16))
e2cS=np.zeros((16,16))
e2cN=np.zeros((16,16))

#b1
#pol_tot=-0.016870524
#b2
pol_tot=-0.021123435

f=open("temp1",'r')
lines1=f.readlines()
f.close()

f=open("temp2",'r')

e1=np.zeros((16,))
e1_0=np.zeros((16,))
e1c=np.zeros((16,))
err1=np.zeros((16,))
e2=np.zeros((16,16))
e2_0=np.zeros((16,16))
e2c=np.zeros((16,16))
e2cS=np.zeros((16,16))
e2cN=np.zeros((16,16))

#b1
#pol_tot=-0.016870524
#b2
pol_tot=-0.021123435

f=open("temp1",'r')
lines1=f.readlines()
f.close()

f=open("temp2",'r')
lines2=f.readlines()
f.close()

f=open("temp",'r')
lines=f.readlines()
f.close()

for line in lines:
   i,j,n,t1,t2,t3,t4,c2,efp,t5,t6=line.split()
   i,j,n,c2,pol=int(i),int(j),n[0],float(c2),float(efp)
   i,j=i-1,j-1
   if n=='N': e2cN[i,j]=pol
   if n=='S': e2cS[i,j]=c2

e1_sum=0.
e1_0c_sum=0.
err_sum=0.
for line1 in lines1:
   fname=line1[:-1]+"_en"
   f=open(fname,'r')
   lines0=f.readlines()
   f.close()
   t,i,j,k=line1.split("_")
   i,j,k=map(int,(i,j,k))
   i=i-1
   e1[i]=float(lines0[0].split()[1])
   e1_0[i]=float(lines0[0].split()[0])
   err1[i]=float(lines0[0].split()[2])
   e1_sum+=e1[i]
   e1_0c_sum+=e1[i]-e1_0[i]
   err_sum+=err1[i]

e2_0c_sum=0.
for line2 in lines2:
   fname=line2[:-1]+"_en"
   f=open(fname,'r')
   lines0=f.readlines()
   f.close()
   t,i,j,k=line2.split("_")
   i,j,k=map(int,(i,j,k))
   i,j=i-1,j-1
   e2[i,j]=float(lines0[0].split()[1])
   e2_0=float(lines0[0].split()[0])
   err2=float(lines0[0].split()[2])
   #e2c[i,j]=e2[i,j]-e1[i]-e1[j]+e2cN[i,j]
   e2c[i,j]=e2[i,j]-e1[i]-e1[j]
   e2_0c_sum+=e2c[i,j]-(e2_0-e1_0[i]-e1_0[j])
   err_sum+=err2+err1[i]+err1[j]

n=16

e2c_sum=0.0
for i in range(1,n):
   for j in range(i):
      #e2c_sum+=e2c[i,j]+e2cS[i,j]
      e2c_sum+=e2c[i,j]

#e,err=e1_sum+e2c_sum+pol_tot,err_sum
e,err=e0eH+e1_0c_sum+e2_0c_sum,err_sum
print "%.2f +/- %.2f"%(round(H2kcal_mol*e,2),round(H2kcal_mol*err,2))
print "%.2f +/- %.2f"%(round(H2kcal_mol*(e-e0H),2),round(H2kcal_mol*err,2))
print "%.2f  %.2f"%(round(H2kcal_mol*e0H,2),round(H2kcal_mol*e0eH,2))
##print H2kcal_mol*(e-e0),H2kcal_mol*(err+err0)
##print e1_0c_sum,e2_0c_sum,e1_0c_sum+e2_0c_sum
#print "%.2f +/- %.2f"%(round(H2kcal_mol*e,2),round(H2kcal_mol*err,2))
#print "%.2f +/- %.2f"%(round(H2kcal_mol*(e0H-e),2),round(H2kcal_mol*err,2))
#print "%.2f +/- %.2f"%(round(H2kcal_mol*e0V,2),round(H2kcal_mol*err0V,2))
#print "%.2f +/- %.2f"%(round(H2kcal_mol*(e0H-e0V),2),round(H2kcal_mol*err0V,2))
#print "%.2f +/- %.2f"%(round(H2kcal_mol*e0,2),round(H2kcal_mol*err0,2))
#print "%.2f +/- %.2f"%(round(H2kcal_mol*(e0H-e0),2),round(H2kcal_mol*err0,2))
