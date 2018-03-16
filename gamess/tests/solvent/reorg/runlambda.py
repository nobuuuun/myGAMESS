#!/usr/bin/python
#
# March, 2016
#
# @author  Soumya Ghosh
#
# @Usage
#
#  runlambda.py arg1 arg2 arg3 arg4 arg5
#
#  where
#  arg1 = name of the 1st equilibrium inputfile without .inp extension
#  arg2 = name of the 2nd equilibrium inputfile without .inp extension
#  arg3 = name of the 1st non-equilibrium inputfile without .inp 
#         extension. This inputfile will be created by this script provided 
#         the equilibrium calculations have finished properly.
#  arg4 = name of the 2nd non-equilibrium inputfile without .inp 
#         extension. This inputfile will be created by this script provided
#         the equilibrium calculations have finished properly.
#  arg5 = name of the textfile without .txt extension that containes
#         summary of the results
#
#
# @Description 
#
#  This script allows the user to 1) run equilibrium calculations
#  2) check *.log files for convergence 3) generate inputfiles for
#  non-equilibrium calculations 4) run the non-equilibrium calculations 
#  5) check the output files for convergence 6) combine the equilibrium
#  and non-equilibrium free energies to print out the final solvent
#  reorganization energy
#
#  USER SPECIFIC NOTES
#
#  This script allows an user to run jobs interactively by directly executing 
#  the 'rungms' script provided with GAMESS. Note that location of 'rungms' 
#  will be different for different users. Moreover, one can manipulate this
#  script to run a submission script, which in turn will call 'rungms'.
#
#
import sys
import os
import subprocess

eql1 = sys.argv[1]
eql2 = sys.argv[2]
neql1 = sys.argv[3]
neql2 = sys.argv[4]
res_reorg = sys.argv[5]

# run checks to confirm that the inputfiles are consistent with each other
# GAMESS groups to check:(1) $DFT (2) $SCF (3) $BASIS (4) $PCM 
# (5) $TESCAV (6) $PCMCAV (7) $REORG (8) $CONTRL (needs special care)
# (9) $DATA (needs special care)

def check_inpfiles(file1,file2,group):

   import difflib

   file1 = file1 + ".inp"
   file2 = file2 + ".inp"
   fileo1 = group[1:] + "_temp1.txt"
   fileo2 = group[1:] + "_temp2.txt"

   if group == "$DATA":
      with open(fileo1,"w") as outfile:
         with open(file1,"r") as infile:
            for line in infile:
               if group in line:
                  for line in infile:
                     if "$END" in line:
                        break
                     else:
                        outfile.write(line.lstrip())

      with open(fileo1,"r") as fin:
         data = fin.read().splitlines(True)
      with open(fileo1,"w") as fout:
         fout.writelines(data[2:]) 
   else: 
      with open(fileo1,"w") as outfile:
         with open(file1,"r") as infile:
            for line in infile:
               if group in line:
                  line1 = line.strip()
                  if "$END" in line1:
                     splitline1 = line1.split()
                     for item in splitline1:
                        if group in item:
                           continue
                        elif "$END" in item:
                           continue
                        else:
                           outfile.write("%s\n" %item)
                     break
                  else:
                     for line in infile:
                        if "$END" in line:
                           break
                        else:
                           splitline1 = line.strip().split()
                           for item in splitline1:
                              outfile.write("%s\n" %item)

   if group == "$CONTRL":
      os.rename(fileo1,"temp1_temp.txt")
      with open(fileo1,"w") as outfile:
         with open("temp1_temp.txt","r") as infile:
            for line in infile:
               if "ICHARG" in line:
                  continue
               elif "MULT" in line:
                  continue
               elif "SCFTYP" in line:
                  continue
               else:
                  outfile.write(line)
      os.remove("temp1_temp.txt")

   if group == "$DATA":
      with open(fileo2,"w") as outfile:
         with open(file2,"r") as infile:
            for line in infile:
               if group in line:
                  for line in infile:
                     if "$END" in line:
                        break
                     else:
                        outfile.write(line.lstrip())

      with open(fileo2,"r") as fin:
         data = fin.read().splitlines(True)
      with open(fileo2,"w") as fout:
         fout.writelines(data[2:]) 
   else: 
      with open(fileo2,"w") as outfile:
         with open(file2,"r") as infile:
            for line in infile:
               if group in line:
                  line1 = line.strip()
                  if "$END" in line1:
                     splitline1 = line1.split()
                     for item in splitline1:
                        if group in item:
                           continue
                        elif "$END" in item:
                           continue
                        else:
                           outfile.write("%s\n" %item)
                     break
                  else:
                     for line in infile:
                        if "$END" in line:
                           break
                        else:
                           splitline1 = line.strip().split()
                           for item in splitline1:
                              outfile.write("%s\n" %item)

   if group == "$CONTRL":
      os.rename(fileo2,"temp2_temp.txt")
      with open(fileo2,"w") as outfile:
         with open("temp2_temp.txt","r") as infile:
            for line in infile:
               if "ICHARG" in line:
                  continue
               elif "MULT" in line:
                  continue
               elif "SCFTYP" in line:
                  continue
               else:
                  outfile.write(line)
      os.remove("temp2_temp.txt")

   filediff = group[1:] + "_difffile.txt"
   with open(filediff,"w") as outfile:
      with open(fileo1,"r") as fileog1:
         with open(fileo2,"r") as fileog2:   
            diff = difflib.unified_diff(fileog1.readlines(),fileog2.readlines(),fromfile = fileog1,tofile = fileog2,n=0,lineterm="")
            for line in diff:
               for prefix in ("---","+++","@@"):
                  if line.startswith(prefix):
                     break
               else:
                  outfile.write(line)
   
   filediffm = [line[1:] for line in open(filediff,"r").readlines() if line[0] == "-"]
   filediffp = [line[1:] for line in open(filediff,"r").readlines() if line[0] == "+"]
   with open(filediff,"w") as fout:
      if len(filediffm) > len(filediffp):
         for line in filediffm:
            if line not in filediffp:
               fout.writelines(line)
      elif len(filediffm) < len(filediffp): 
         for line in filediffp:
            if line not in filediffm:
               fout.writelines(line)
      elif len(filediffm) == len(filediffp):
         for line in filediffm:
            if line not in filediffp:
               fout.writelines(line)
         for line in filediffp:
            if line not in filediffm:
               fout.writelines(line)
   filedifflist = [line for line in open(filediff,"r").readlines()] 
   os.remove(fileo1)      
   os.remove(fileo2)      
   os.remove(filediff)
   if len(filedifflist) == 0:
#      print ("True_" + group[1:])
      return "True"
   else:
#      print ("False_" + group[1:])
      return "False"

   
chk_inp_contrl = check_inpfiles(eql1,eql2,"$CONTRL")
chk_inp_scf = check_inpfiles(eql1,eql2,"$SCF")
chk_inp_dft = check_inpfiles(eql1,eql2,"$DFT")
chk_inp_basis = check_inpfiles(eql1,eql2,"$BASIS")
chk_inp_pcm = check_inpfiles(eql1,eql2,"$PCM")
chk_inp_pcmcav = check_inpfiles(eql1,eql2,"$PCMCAV")
chk_inp_tescav = check_inpfiles(eql1,eql2,"$TESCAV")
chk_inp_reorg = check_inpfiles(eql1,eql2,"$REORG")
chk_inp_data = check_inpfiles(eql1,eql2,"$DATA")

if (chk_inp_contrl == "True" and 
    chk_inp_scf == "True" and 
    chk_inp_dft == "True" and 
    chk_inp_basis == "True" and
    chk_inp_pcm == "True" and
    chk_inp_pcmcav == "True" and
    chk_inp_tescav == "True" and
    chk_inp_reorg == "True" and
    chk_inp_data == "True"):
   print ("All checks complete. The script will now submit the equilibrium jobs.\n Modify your rungms file to make sure that the *.log and *.dat files are\n in the same directory as this script. Otherwise modify this python script accordingly.")  
# submit eq jobs
   subeqfile1 = ("/home/ghoshs/gamess/gms_dev/GMS_ORIGIN/GAMESS_REORG/rungms " + eql1 + " 00 1 >& " + eql1 + ".log")
   subprocess.Popen(subeqfile1,shell=True,stdout=None)
   subeqfile2 = ("/home/ghoshs/gamess/gms_dev/GMS_ORIGIN/GAMESS_REORG/rungms " + eql2 + " 00 1 >& " + eql2 + ".log")
   subprocess.Popen(subeqfile2,shell=True,stdout=None)
else:
   print ("Check inputfiles for keyword inconsistency")

fileq1log = eql1 + ".log"
fileq1dat = eql1 + ".dat"
fileq2log = eql2 + ".log"
fileq2dat = eql2 + ".dat"

# run checks to confirm that the jobs are finished and converged

def check_subfiles(file1):
   
   import time

   time.sleep(2*60)

   i=1
   while i == 1:
      fileo1 = open(file1).readlines()        
      for line in reversed(fileo1):
         if "EXECUTION OF GAMESS TERMINATED NORMALLY" in line:
            f1stat = "Calculation finished"
            break
         elif "EXECUTION OF GAMESS TERMINATED -ABNORMALLY-" in line:
            f1stat = "Calculation died"
            break
         else:
            f1stat = "Calculation still running"         

      if (f1stat == "Calculation finished"):
         for line in reversed(fileo1):
            if "DENSITY CONVERGED" in line:
               file1stat = "Calculation finished and converged"
               break
            else:
               file1stat = "Calculation finished but not converged"
         return file1stat
         break
      elif (f1stat == "Calculation died"):
         file1stat = f1stat
         return file1stat
         break
   

eq1logstat = check_subfiles(fileq1log)
print "Equilibrium1 " + eq1logstat
eq2logstat = check_subfiles(fileq2log)
print "Equilibrium2 " + eq2logstat

def generate_neqfiles(file1,file2,file3):
   import shutil

   file1 = file1 + ".inp"
   file3 = file3 + ".inp"
   shutil.copy2(file1,file3)
   fileo1 = file1[:-4] + ".dat"
   fileo2 = file2 + ".dat"
   
   fneq = open(file3,'a')
   fo1 = open(fileo1,'r')
   fo2 = open(fileo2,'r')
   with fneq as outfile:
      with fo1 as infile:
         for line in infile:
            if "$SRFCHG" in line:
               outfile.write(line)
               for line in infile:
                  if "$END" in line:
                     break
                  else:
                     outfile.write(line)
      fo1.close()
      outfile.write(" $END\n")
      with fo2 as infile:
         for line in infile:
            if "$SRFCHG" in line:
               outfile.write(line)
               for line in infile:
                  if "$END" in line:
                     break
                  else:
                     outfile.write(line)
      fo2.close()
      outfile.write(" $END")
   fneq.close()        
   fneq = open(file3,'r')
   for line in fneq:
      if "$REORG" in line:
         strreorg = line.strip()
         if "$END" in strreorg:
            strtrim = strreorg[6:-4]
            strnew = "$REORG\n" + strtrim + "\n RLMIT=SC\n $END"
         else:
            strnew = strreorg + "\n RLMIT=SC"
   fneq.close()
   fneq = open(file3).read()
   fneq = fneq.replace(strreorg,strnew)   
   fneq = fneq.replace("IPCHG","IRCHG")
   fneqn = open(file3,'w')
   fneqn.write(fneq)
   fneqn.flush()
   fneqn.close()

# generate non-equilibrium files by combining equilibrium ("ox" and "red") .dat and .inp files
 
if (eq1logstat == "Calculation finished and converged") and (eq2logstat == "Calculation finished and converged"):
   generate_neqfiles(eql1,eql2,neql1)
   generate_neqfiles(eql2,eql1,neql2)

# submit neq jobs
   
   subneqf1 = ("/home/ghoshs/gamess/gms_dev/GMS_ORIGIN/GAMESS_REORG/rungms " + neql1 + " 00 1 >& " + neql1 + ".log")
   subprocess.Popen(subneqf1,shell=True,stdout=None)
   subneqf2 = ("/home/ghoshs/gamess/gms_dev/GMS_ORIGIN/GAMESS_REORG/rungms " + neql2 + " 00 1 >& " + neql2 + ".log")
   subprocess.Popen(subneqf2,shell=True,stdout=None)

# check noneq logfiles

   fileneq1log = neql1 + ".log"
   fileneq2log = neql2 + ".log"
   neq1logstat = check_subfiles(fileneq1log)
   print "Non-equilibrium1 "+neq1logstat
   neq2logstat = check_subfiles(fileneq2log)
   print "Non-equilibrium2 "+neq2logstat
else:
   print "Check logfiles"

# calculate solvent reorganization energy and store it in results file
 
if (neq1logstat == "Calculation finished and converged") and (neq2logstat == "Calculation finished and converged"):
   res_reorg_txt = res_reorg + ".txt"
   fresults = open(res_reorg_txt,"w")
   feqinp1 = eql1+".inp"
   feqinpo1 = open(feqinp1).readlines()
   feqinp2 = eql2+".inp"
   feqinpo2 = open(feqinp2).readlines()
   fneqo1 = open(fileneq1log).readlines()
   feqo1 = open(fileq1log).readlines()
   fneqo2 = open(fileneq2log).readlines()
   feqo2 = open(fileq2log).readlines()
   with open("temp1.txt","w") as outfile:
      with open(feqinp1,"r") as infile:
         for line in infile:
            if "$DATA" in line:
               for line in infile:
                  if "$END" in line:
                     break
                  else:
                     outfile.write(line)
#
   with open("temp1.txt","r") as fin:
      data = fin.read().splitlines(True)
   with open("temp1.txt","w") as fout:
      fout.writelines(data[2:]) 

   with fresults as outfile:
      outfile.write("\n          __                                 ")
      outfile.write("\n         /  \                                ")
      outfile.write("\n         \__   __             __      ___    ")
      outfile.write("\n            \ /  \ |   \   / |__ |\ |  |     ")
      outfile.write("\n         \__/ \__/ |__  \_/  |__ | \|  |     \n")
      outfile.write("\n               ___                            ") 
      outfile.write("\n              |   |                           ")
      outfile.write("\n              |_ _|   __  __   _   __    _          __   _   ___    __       ")
      outfile.write("\n              |   \  |__ /  \ |_] / __  /_\  |\ | |  /  /_\   |  | /  \ |\ | ")
      outfile.write("\n              |    \ |__ \__/ | \ \__/ /   \ | \| | /_ /   \  |  | \__/ | \| \n")
      outfile.write("\n                     ___                                                     ")
      outfile.write("\n                    |                                                        ")
      outfile.write("\n                    |___       __  _   __                                    ")
      outfile.write("\n                    |    |\ | |__ |_] / __ \_/                               ")
      outfile.write("\n                    |___ | \| |__ | \ \__/  |                                \n\n")
      outfile.write("...............................................................................")
      outfile.write("\n\n This script prints out the final average solvent reorganization energy at a ")
      outfile.write("\n particular solute geometry as implemented within the framework of IEF-PCM in")
      outfile.write("\n GAMESS quantum chemistry software package.                                  ")
      outfile.write("\n\n The solvent reorganization energy is an important parameter in Marcus     ")
      outfile.write("\n Theory that connects the free energy of activation and the equilibrium    ")
      outfile.write("\n reaction free energy for electron transfer and proton-coupled electron    ")
      outfile.write("\n transfer reactions. This parameter is a measure of the energy penalty      ")
      outfile.write("\n required to change the equilibrium solvent configuration upon electron    ")
      outfile.write("\n transfer.                                                                   \n\n")
      outfile.write("...............................................................................")
      outfile.write("\n\n                  REFERENCES                 \n")
      outfile.write("\n    1. Ghosh,S.; Horvath,S.; Soudackov,A.V.; Hammes-Schiffer,S.              ")
      outfile.write("\n       Electrochemical Solvent Reorganization Energies in the Framework of   ")                                           
      outfile.write("\n       the Polarizable Continuum Model                                       ")
      outfile.write("\n       J. Chem. Theory Comput. 2014, 10, 2091-2102.                          \n")
      outfile.write("\n    2. Ghosh,S.; Hammes-Schiffer,S.                                          ")
      outfile.write("\n       Calculation of Electrochemical Reorganization Energies for Redox      ")
      outfile.write("\n       Molecules at Self-Assembled Monolayer Modified Electrodes             ")
      outfile.write("\n       J. Phys. Chem. Lett. 2015, 6, 1-5.                                    \n\n")
      outfile.write("...............................................................................")
      outfile.write("\n               ================              ")
      outfile.write("\n                INPUT KEYWORDS               ")
      outfile.write("\n               ================              \n")
      outfile.write("\n     ------------- STATE 1 -------------     \n")
      for line in feqinpo1:
         if "$DATA" in line:
            break
         else:
            outfile.write(line)
      outfile.write("\n     ------------- STATE 2 -------------     \n")
      for line in feqinpo2:
         if "$DATA" in line:
            break
         else:
            outfile.write(line)
      outfile.write("\n               =================              ")
      outfile.write("\n                XYZ-COORDINATES               ")
      outfile.write("\n               =================              \n\n")
      with open("temp1.txt","r") as infile:
         for line in infile:
            outfile.write(line)
      os.remove("temp1.txt")
      outfile.write("...............................................................................")
      outfile.write("\n                                       |                                       ")
      outfile.write("\n                                   =========                                   ")
      outfile.write("\n                                   |RESULTS|                                   ")
      outfile.write("\n                                   =========                                   \n\n")
      with open(fileq1log,"r") as infile:
         for line in infile:
            if "INPUT FOR PCM SOLVATION CALCULATION" in line:
               for line in infile:
                  if "EPS" in line:
                     line_eps = line.lstrip()
                     eps_static = float(line_eps[27:35])
                     eps_inf = float(line_eps[47:55])
                     eps_data = "        (SOLVENT PROPERTIES:: EPS_STATIC = " + str(eps_static) + "; EPS_INFINITY = " + str(eps_inf) + ")\n"
                     outfile.write(eps_data)
                     break
               break
      outfile.write("\n                      ------------- STATE 1 -------------                    \n")
      for line in reversed(fneqo1):
         if "NONEQ FREE ENG IN SOLVENT  =   <PSI|H(0)+V/2|PSI>   =" in line:
            fneq1str = line.strip()
            fneq1flt = float(fneq1str[-25:-6])
            fneq1res = "\n NONEQUILIBRIUM FREE ENERGY IN SOLVENT  =   " + str(fneq1flt) + " A.U.\n"
            outfile.write(fneq1res)
      for line in reversed(feqo1):
         if "FREE ENERGY IN SOLVENT = <PSI| H(0)+V/2 |PSI>       =" in line:
            feq1str = line.strip()
            feq1flt = float(feq1str[-25:-6])
            feq1res = " EQUILIBRIUM FREE ENERGY IN SOLVENT     =   " + str(feq1flt) + " A.U.\n"
            outfile.write(feq1res)
      outfile.write("                                          ------------------------\n")
      out_reorg_eng1 = (fneq1flt - feq1flt)*27.211396
      lambda1_res = " THE SOLVENT REORGANIZATION ENERGY      =   " + "{0:.6f}".format(out_reorg_eng1) + " eV\n"
      outfile.write(lambda1_res)
      outfile.write("\n                      ------------- STATE 2 -------------                    \n")
      for line in reversed(fneqo2):
         if "NONEQ FREE ENG IN SOLVENT  =   <PSI|H(0)+V/2|PSI>   =" in line:
            fneq2str = line.strip()
            fneq2flt = float(fneq2str[-25:-6])
            fneq2res = "\n NONEQUILIBRIUM FREE ENERGY IN SOLVENT  =   " + str(fneq2flt) + " A.U.\n"
            outfile.write(fneq2res)
      for line in reversed(feqo2):
         if "FREE ENERGY IN SOLVENT = <PSI| H(0)+V/2 |PSI>       =" in line:
            feq2str = line.strip()
            feq2flt = float(feq2str[-25:-6])
            feq2res = " EQUILIBRIUM FREE ENERGY IN SOLVENT     =   " + str(feq2flt) + " A.U.\n"
            outfile.write(feq2res)
      outfile.write("                                          ------------------------\n")
      out_reorg_eng2 = (fneq2flt - feq2flt)*27.211396
      lambda2_res = " THE SOLVENT REORGANIZATION ENERGY      =   " + "{0:.6f}".format(out_reorg_eng2) + " eV\n"
      outfile.write(lambda2_res)
      outfile.write("\n*******************************************************************************")
      out_reorg_eng = (out_reorg_eng1 + out_reorg_eng2)/2.0
      lambda_res = "\n THE AVERAGE SOLVENT REORGANIZATION ENERGY FOR THIS STRUCTURE  =   " + "{0:.6f}".format(out_reorg_eng) + " eV\n\n"
      outfile.write(lambda_res)
   fresults.close()   
else:
   print "Check non-eq logfiles"

