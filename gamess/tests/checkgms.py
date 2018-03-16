#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import json
import sys

from checkgms_utils import *
from checkgms_stable import *

"""
If --json_create flag is passed to trigger the creation of new *.json validation files then it will look or log
files containg the extension specified by file_extension to parse in order to generate the new *.json validation file.
"""
file_extension=".VALIDATION"

"""
In some cases we may have a *.json validation file that we wish to protect from being over-written even when the --json_create
flag is passed in. A way to allow this is by renaming the file so that it ends with what is specified by json_protect.

Meaning if --json_protect is passed in then the script will skip all input files that match the filename patter *project.json
"""
json_protect="protect.json"

run_arguments=parse_arguments()

print(c_box("Run parameters"))
print (json.dumps(run_arguments,indent=2))

#Populate the parser_groups array
parse_groups=[]
parse_groups=get_parse_groups(run_arguments=run_arguments)

#Populate the log_file_paths array
log_file_paths=[]
if run_arguments["json_create"]:
  log_file_paths=get_log_file_paths(folder_string_match=run_arguments["filter_folder"],file_string_match=run_arguments["filter_file"],folder_string_skip=run_arguments["skip_folder"],file_string_skip=run_arguments["skip_file"],file_extension=file_extension)
else:
  log_file_paths=get_log_file_paths(folder_string_match=run_arguments["filter_folder"],file_string_match=run_arguments["filter_file"],folder_string_skip=run_arguments["skip_folder"],file_string_skip=run_arguments["skip_file"])

#Loop through the log_file_paths array and validate
for filenum, log_file_path in enumerate(log_file_paths,start=1):

  if run_arguments["json_create"]:
    validation_result=checkgms(filenum=filenum,log_file_path=log_file_path,log_file_count=len(log_file_paths),run_arguments=run_arguments,parse_groups=parse_groups,file_extension=file_extension,json_protect=json_protect)
  else:
    validation_result=checkgms(filenum=filenum,log_file_path=log_file_path,log_file_count=len(log_file_paths),run_arguments=run_arguments,parse_groups=parse_groups)

  if "skip" in str(validation_result):
    continue

  if run_arguments["exit_on_fail"]:
    if "pass" not in str(validation_result):
      print("Exiting on failure!")
      sys.exit(1)
