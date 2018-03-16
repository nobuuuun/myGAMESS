#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import json
import sys
import os
import time
import subprocess
import mmap
import datetime

from checkgms_utils import *
from checkgms_stable import *

codecoverage=[]

total_coverage=0
total_linesofcode=0

first=True

run_arguments=parse_arguments(validation=False)
rungms_path="misc/automation/rungms"

"""
print(c_box("Run parameters"))
print (json.dumps(run_arguments,indent=2))
"""

script_path = os.path.dirname(os.path.realpath(__file__))

#Populate the log_file_paths array
input_file_paths=[]
input_file_paths=get_input_file_paths(folder_string_match=run_arguments["filter_folder"],file_string_match=run_arguments["filter_file"],folder_string_skip=run_arguments["skip_folder"],file_string_skip=run_arguments["skip_file"],script_path=script_path)

#Loop through the log_file_paths array and validate
for filenum, input_file_path in enumerate(input_file_paths,start=1):

  short_input_file_path=input_file_path.split("/tests/",1)[-1]

  #Search for prescence of "TRAVIS-CI SKIP"
  with open(input_file_path,'r',encoding="utf-8",errors='ignore') as opened_file:
    parse_memory_map = mmap.mmap(opened_file.fileno(), 0, access=mmap.ACCESS_READ)
    regex_string="TRAVIS-CI SKIP"
    regex=re.compile(str.encode(regex_string,'ascii'), re.MULTILINE)
    match = regex.search(parse_memory_map)
    #If found then skip
    if match:
      if run_arguments["debug"]:
        print(l_box_small("Skipping input file"),file_progress(filenum,len(input_file_paths)),short_input_file_path)
      continue

    #If --test_type is passed in:
    if "small" in run_arguments["test_type"]:
      regex_string="TRAVIS-CI SMALL"
    elif "medium" in run_arguments["test_type"]:
      regex_string="TRAVIS-CI MEDIUM"
    elif "large" in run_arguments["test_type"]:
      regex_string="TRAVIS-CI LARGE"

    if len(run_arguments["test_type"]) > 0:
      regex=re.compile(str.encode(regex_string,'ascii'), re.MULTILINE)
      match = regex.search(parse_memory_map)
      if not match:
        if run_arguments["debug"]:
          print(l_box_small("Skipping input file"),file_progress(filenum,len(input_file_paths)),short_input_file_path)
        continue

  try:

    run_command=rungms_path+" "+input_file_path+" 00 "+run_arguments["ncpus"]+" "+run_arguments["ncpus"]+ " > "+input_file_path.replace(".inp",".log")+" 2>&1"

    if not run_arguments["coverage"]:
        print(c_box_small(datetime.datetime.now().time()),l_box_small("Running input file"),file_progress(filenum,len(input_file_paths)),short_input_file_path)

    if run_arguments["stderr"]:
      sys.stderr.write(l_box_small("Running input file")+" "+file_progress(filenum,len(input_file_paths))+" "+short_input_file_path+"\n")
      sys.stderr.flush()

    if run_arguments["debug"] or run_arguments["dryrun"]:
      print(run_command)

    if not run_arguments["dryrun"]:
      os.system(run_command)

    if run_arguments["coverage"]:
      result = subprocess.run(["tests/coverage"], stdout=subprocess.PIPE).stdout.decode('utf-8')

      lines=result.splitlines()

      if first:
        print(c_large("Test Input"),c_small("% Total Coverage"))
        print(seperator(columns=83))
        codecoverage.append(["Source Files"])
        codecoverage.append(["Lines of Code"])

      codecoverage.append([short_input_file_path])

      for line in range(len(lines)):
        if "File" in lines[line]:
          sourcefile=lines[line].split("'")[-2]
          coverage=lines[line+1].split(":")[-1].split("%")[0]
          linesofcode=lines[line+1].split(" ")[-1]

          if "No executable lines" in lines[line+1]:
            coverage="0.00"
            linesofcode=0

          if first:
            codecoverage[0].append(sourcefile)
            codecoverage[1].append(linesofcode)

          codecoverage[-1].append(coverage)

      if first:
        codecoverage[0].append("Total")
        loc_sum=0
        for loc in codecoverage[1][1:-1]:
          loc_sum=loc_sum+int(loc)
        codecoverage[1].append(loc_sum)

      codecoverage[-1].append(lines[-1].split(":")[-1].split("%")[0])

      if first:
        first=False

      if float(codecoverage[-1][-1]) > float(codecoverage[-2][-1]):
        print(l_box_large(codecoverage[-1][0]),r_small(codecoverage[-1][-1]))
      else:
        if float(codecoverage[-2][-1]) < 100.0:
          print(l_box_large(codecoverage[-1][0]),r_small(codecoverage[-1][-1]),c_small("X"))
        else:
          print(l_box_large(codecoverage[-1][0]),r_small(codecoverage[-1][-1]))

    #To permit a keyboard interrupt
    time.sleep(2)

  except KeyboardInterrupt:
    sys.exit(1)
  except:
    sys.exit(1)

if run_arguments["coverage"]:
  with open("tests/coverage.dat","w") as coverage_file:
    for row in codecoverage:
     coverage_file.write(",".join([ str(x) for x in row]))
     coverage_file.write("\n")
     pass

  print(seperator(columns=83))
  print(l_box_large("Total % Coverage"),r_small(codecoverage[-1][-1]))
  print(l_box_large("Total Lines of Code"),r_small(codecoverage[1][-1]))
