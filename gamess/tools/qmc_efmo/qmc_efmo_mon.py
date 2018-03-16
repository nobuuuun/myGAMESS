#!/usr/bin/env python

from fireworks import Firework, LaunchPad, ScriptTask, FWorker
from fireworks.queue.queue_launcher import launch_rocket_to_queue, rapidfire
from fireworks.utilities.fw_serializers import load_object_from_file

import os

launchpad = LaunchPad(host="mongodb03", port=27017, name="fw_gamess_qmcpack", username=<password>, password=<username>)
launchpad.reset('2017-09-14')

f=open("temp1",'r')
lines=f.readlines()
f.close()

n=len(lines)

for line in lines:
   frag=line[:-1]
   firetask=ScriptTask.from_str('cp ../../../'+frag+'.inp .;'+\
   'srun -n 1 /global/homes/f/fzahari/Work/Programs/gamess_edison/rungms '+frag+' > gamess.out;'+\
   'srun -n 1 /global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/build_edison/bin/convert4qmc -gamessAscii gamess.out;'+\
   'cp sample.Gaussian-G2.xml h2o.wfs1.xml;'+\
   'cp sample.Gaussian-G2.ptcl.xml h2o.ptcl.xml;'+\
   'cp ../../../templ_files/h2o.optm1.xml .;'+\
   'cp ../../../templ_files/O.xml .;'+\
   'cp ../../../templ_files/H.xml .;'+\
   'export OMP_NUM_THREADS=1;'+\
   'srun -n 528 /global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/build_edison/bin/qmcpack h2o.optm1.xml &> h2o.optm1.out;'+\
   '../../../scripts/jopt1.py;'+\
   '../../../scripts/j3.py;'+\
   'cp ../../../templ_files/h2o.optm2.xml .;'+\
   'srun -n 528 /global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/build_edison/bin/qmcpack h2o.optm2.xml &> h2o.optm2.out;'+\
   '../../../scripts/jopt2.py;'+\
   'cp ../../../templ_files/h2o.dmc.xml .;'+\
   'srun -n 528 /global/homes/f/fzahari/Work/Programs/qmcpack-3.0.0/build_edison/bin/qmcpack h2o.dmc.xml &> h2o.dmc.out;'+\
   '../../../scripts/en.py '+frag+' > ../../../'+frag+'_en;')
   firework = Firework(firetask)
   launchpad.add_wf(firework)

if not os.path.exists('mon'):
   os.makedirs('mon')

template_file="/global/homes/f/fzahari/queue_tests/my_qadapter.yaml"
qadapter=load_object_from_file(template_file)
rapidfire(launchpad,FWorker(),qadapter,launch_dir="mon",njobs_queue=n)
