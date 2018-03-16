import json
import os

import mmap

from checkgms_utils import *
from checkgms_parsers import parser

def checkgms(filenum=None,log_file_path=None,log_file_count=0,run_arguments={},parse_groups=None,file_extension=".log",json_protect=None):

  if log_file_path is None:
    print(lr_box("log_file_path is undefined!"))
    return

  if parse_groups is None:
    print(lr_box("parse_groups is undefined!"))
    return

  parsed_json={}
  validation_json={}
  validated_json={}

  log_filename=log_file_path.split("/")[-1]
  base_filename=log_filename.replace(file_extension,"")
  validation_file_path=log_file_path.replace(file_extension,".json")

  with open(log_file_path,'r',encoding="utf-8",errors='ignore') as opened_file:
    print(l_box("Parsing"),file_progress(filenum,log_file_count),l_(log_file_path))
    if not run_arguments["dryrun"]:
      #Create a memory map of the log file for faster searching if a given regex_string exists within the document
      parse_memory_map = mmap.mmap(opened_file.fileno(), 0, access=mmap.ACCESS_READ)
      #Parse the log file and return a JSON representing the parsed log file
      values_json=parse(file_handle=opened_file,parse_memory_map=parse_memory_map,run_arguments=run_arguments,parse_groups=parse_groups)
    else:
      values_json={}

    parsed_json["name"]=log_filename.replace(file_extension,".log")
    parsed_json["validation"]=values_json

    if run_arguments["verbose_parsing"]:
      #Formatted printing of the parsed JSON
      print_parsed_JSON(parsed_json=parsed_json,run_arguments=run_arguments)

    #Check if validation file exists. Yes, validate. No, save the parsed JSON as a validation file.
    if os.path.isfile(validation_file_path) and not run_arguments["json_create"]:
      print(l_box("Validating"),file_progress(filenum,log_file_count),l_(log_file_path))
      if not run_arguments["dryrun"]:
        with open(validation_file_path,'r',encoding="utf-8",errors='ignore') as validation_file:
          validation_json=json.load(validation_file)
        #Validate parsed JSON and return a modified JSON containg calculated errors and validation results ("pass/fail/skip/validation")
        validated_json=validate(validation_json=validation_json,parsed_json=parsed_json)
        if run_arguments["verbose_validation"] and validated_json["result"]!="validation":
          #Formatted printing of the validated JSON
          print_validated_JSON(validated_json=validated_json,run_arguments=run_arguments)
      else:
        #Early return if its a dry run
        if run_arguments["dryrun"]:
          validated_json["result"]="skip"
          return validated_json["result"]
      #Simplified formatted printing for validation results
      if validated_json["result"] == "pass":
        print(l_box("Validation result"),file_progress(filenum,log_file_count),l_(log_file_path),pass_box())
      else:
        print(l_box("Validation result"),file_progress(filenum,log_file_count),l_(log_file_path),fail_box())
      return validated_json["result"]
    #Check if validation file exists AND if it is projected
    elif os.path.isfile(validation_file_path) and json_protect is not None and json_protect in validation_file_path:
      validated_json["result"]="skip"
      print(l_box("Skipping validation file create"),file_progress(filenum,log_file_count),l_(validation_file_path))
      return validated_json["result"]
    #Chekf if we are skipping validation file creation
    elif run_arguments["skip_json_create"]:
      validated_json["result"]="skip"
      print(l_box("Skipping validation file create"),file_progress(filenum,log_file_count),l_(validation_file_path))
      return validated_json["result"]
    else:
      validated_json["result"]="skip"
      print(l_box("Creating validation file"),file_progress(filenum,log_file_count),l_(validation_file_path))
      if not run_arguments["dryrun"]:
        with open(validation_file_path,'w',encoding="utf-8",errors='ignore') as validation_file:
          validation_file.write(json.dumps(parsed_json,indent=2))
      return validated_json["result"]

def validate(validation_json=None,parsed_json=None):
  if validation_json is None:
    print(lr_box("validation_json is undefined!"))
    return

  if parsed_json is None:
    print(lr_box("parsed_json is undefined!"))
    return

  #If filename matches, proceed
  if validation_json["name"] == parsed_json["name"]:
    #If the number of validation values are the same, proceed
    if len(validation_json["validation"]) == len(parsed_json["validation"]):
      parsed_json["result"]="pass"
      #Loop through validation values
      for i in range(len(parsed_json["validation"])):
        #If validation value name matches, proceed
        if validation_json["validation"][i]["name"] == parsed_json["validation"][i]["name"]:

          #Changed: using absolute values during error calculation. Allowing sign changes.
          if hasattr(validation_json["validation"][i]["value"], '__len__') and (not isinstance(validation_json["validation"][i]["value"], str)):
            error=[]
            for validation_element_value, parsed_element_value in zip(validation_json["validation"][i]["value"],parsed_json["validation"][i]["value"]):
              element_error=abs(float(validation_element_value))-abs(float(parsed_element_value))
              error.append(element_error)
          else:
            if validation_json["validation"][i]["type"] == "float":
              error=abs(float(validation_json["validation"][i]["value"]))-abs(float(parsed_json["validation"][i]["value"]))
            elif validation_json["validation"][i]["type"] == "integer":
              error=abs(int(validation_json["validation"][i]["value"]))-abs(int(parsed_json["validation"][i]["value"]))

          parsed_json["validation"][i]["error"]=error

          #Find largest absolute if we are dealing with an array of values and set that as our value for error
          if hasattr(validation_json["validation"][i]["value"], '__len__') and (not isinstance(validation_json["validation"][i]["value"], str)):
            min_error=min(error)
            max_error=max(error)
            if abs(min_error) > abs(max_error):
              error=min_error
            else:
              error=max_error
            parsed_json["validation"][i]["max_error"]=error

          if abs(error) <= validation_json["validation"][i]["tolerance"]:
            parsed_json["validation"][i]["result"]="pass"
          else:
            parsed_json["validation"][i]["result"]="fail"
            parsed_json["result"]="fail"
        else:
          print(lr_box("Validation entries do not match!"),validation_json["validation"][i]["name"],'vs.',parsed_json["validation"][i]["name"])
          parsed_json["validation"][i]["result"]="fail"
          parsed_json["result"]="validation"
    else:
      print(lr_box("Number of validations do not match!"),len(validation_json["validation"]),'vs.',len(parsed_json["validation"]))
      parsed_json["result"]="validation"
  else:
    print(lr_box("Validating incorrect files!"),validation_json["name"],'vs.',parsed_json["name"])
    parsed_json["result"]="validation"

  return parsed_json

def parse(file_handle=None,parse_memory_map=None,run_arguments=None, parse_groups=None):
  parsed_values=[]

  for parse_group_json in parse_groups:
    #Populate parse header or resort to parser group filename
    try:
      parse_header=parse_group_json["parse_header"]
    except:
      parse_header="Parse Group"

    if run_arguments["debug"] : print (lb_box("Debug : ",parse_header))

    #Loop over parse group
    for parse in parse_group_json["parse_group"]:

      #Process parse JSON
      try:
        regex_string=parse["regex_string"]
      except:
        print(lr_box("regex_string was not found!"))
        continue

      try:
        parse_name=parse["parse_name"]
      except:
        parse_name="Value Name"

      try:
        parse_case_sensitive=parse["parse_case_sensitive"]
      except:
        parse_case_sensitive=True

      try:
        parse_type=parse["parse_type"]
      except:
        parse_type="float"

      parse_integer=False
      if parse_type == "integer":
        parse_integer=True

      try:
        parse_precision=parse["parse_precision"]
      except:
        parse_precision=10

      try:
        parse_tolerance=parse["parse_tolerance"]
      except:
        parse_tolerance=None

      try:
        parse_value_index=parse["parse_value_index"]
      except:
        parse_value_index=0

      try:
        parse_line_offset=parse["parse_line_offset"]
      except:
        parse_line_offset=0

      try:
        parse_instance=parse["parse_instance"]
      except:
        parse_instance=-1

      try:
        parse_required=parse["parse_required"]
      except:
        parse_required=False

      try:
        parse_recipe=parse["parse_recipe"]
      except:
        parse_recipe=None

      if parse_recipe is not None:
        if parse_recipe == "lastinstance_firstindex":
          parse_instance=-1
          parse_value_index=0
        elif parse_recipe == "lastinstance_lastindex":
          parse_instance=-1
          parse_value_index=-1
        elif parse_recipe == "firstinstance_firstindex":
          parse_instance=1
          parse_value_index=0
        elif parse_recipe == "firstinstance_lastindex":
          parse_instance=1
          parse_value_index=-1
        else:
          print(lr_box("parse_recipe does not currenty exist!"))
          continue

      #Set parse_multiline to True if we are parsing an array
      parse_multiline=False
      if "array" in parse_type:
        parse_multiline=True

      parse_memory_map_temp=parse_memory_map

      parse_multiple=False
      try:
        parse_multiple_indexes=parse["parse_multiple_indexes"]
        parse_multiple_label=parse["parse_multiple_label"]
        parse_multiple=True
      except:
        parse_multiple_instances=None

      try:
        parse_multiple_instances=parse["parse_multiple_instances"]
        parse_multiple_label=parse["parse_multiple_label"]
        parse_multiple=True
      except:
        parse_multiple_instances=None

      try:
        parse_multiple_lines=parse["parse_multiple_lines"]
        parse_multiple_label=parse["parse_multiple_label"]
        parse_multiple=True
      except:
        parse_multiple_lines=None

      if parse_multiple:

        #Parse multiple instances
        if parse_multiple_instances is not None:
          #Last N instances
          if parse_multiple_instances < 0:
            #Loop over parse instances
            for instance in range(parse_multiple_instances,0):
              parse_name_instance=parse_name+" : "+parse_multiple_label+" "+str((instance))
              parse_instance=instance
              file_handle.seek(0)
              parsed_value=parser(
                regex_string,file_handle,
                parse_memory_map=parse_memory_map_temp,
                parse_value_index=parse_value_index,
                parse_line_offset=parse_line_offset,
                parse_instance=parse_instance,
                parse_integer=parse_integer,
                parse_case_sensitive=parse_case_sensitive,
                parse_multiline=parse_multiline)

              if parsed_value is not None:

                #Calculate precision
                if parse_precision < 0:
                  rhs=parsed_value.split(".")[1]
                  parse_precision=len(rhs)

                #Calculate tolerance
                if parse_tolerance is not None:
                  if parse_tolerance < 0:
                    lhs=parsed_value.replace("-","").split(".")[0]
                    rhs=parsed_value.replace("-","").split(".")[1]
                    if len(lhs) > 4:
                      parse_tolerance=5.0*10**(len(lhs)-12)
                    else:
                      parse_tolerance=5.0*10**(-(len(rhs)-2))

                parsed_values.append(
                  create_JSON(
                    value_name=parse_name_instance,
                    value_value=parsed_value,
                    value_header=parse_header,
                    value_type=parse_type,
                    value_precision=parse_precision,
                    value_tolerance=parse_tolerance))
              else:
                #We don't break because we are looping in reverse
                pass

            if parsed_value is None:
              if parse_required:
                #Break out of the for loop over parse group
                break

          #First N instances
          if parse_multiple_instances > 0:
            for instance in range(1,parse_multiple_instances+1):
              parse_name_instance=parse_name+" : "+parse_multiple_label+" "+str(instance)
              parse_instance=instance
              file_handle.seek(0)
              parsed_value=parser(
                regex_string,file_handle,
                parse_memory_map=parse_memory_map_temp,
                parse_value_index=parse_value_index,
                parse_line_offset=parse_line_offset,
                parse_instance=parse_instance,
                parse_integer=parse_integer,
                parse_case_sensitive=parse_case_sensitive)

              if parsed_value is not None:

                #Calculate precision
                if parse_precision < 0:
                  rhs=parsed_value.split(".")[1]
                  parse_precision=len(rhs)

                #Calculate tolerance
                if parse_tolerance is not None:
                  if parse_tolerance < 0:
                    lhs=parsed_value.replace("-","").split(".")[0]
                    rhs=parsed_value.replace("-","").split(".")[1]
                    if len(lhs) > 4:
                      parse_tolerance=5.0*10**(len(lhs)-12)
                    else:
                      parse_tolerance=5.0*10**(-(len(rhs)-2))

                parsed_values.append(
                  create_JSON(
                    value_name=parse_name_instance,
                    value_value=parsed_value,
                    value_header=parse_header,
                    value_type=parse_type,
                    value_precision=parse_precision,
                    value_tolerance=parse_tolerance))
              else:
                #Break out of the for loop over parse instances
                break

            if parsed_value is None:
              if parse_required:
                #Break out of the for loop over parse group
                break

        elif parse_multiple_lines is not None:
          for line in range(0,parse_multiple_lines):
            parse_name_line=parse_name+" : "+parse_multiple_label+" "+str(line+1)
            parse_line_offset_=parse_line_offset+line
            file_handle.seek(0)
            parsed_value=parser(
              regex_string,file_handle,
              parse_memory_map=parse_memory_map_temp,
              parse_value_index=parse_value_index,
              parse_line_offset=parse_line_offset_,
              parse_instance=parse_instance,
              parse_integer=parse_integer,
              parse_case_sensitive=parse_case_sensitive,
              parse_multiline=parse_multiline)

            if parsed_value is not None:

              #Calculate precision
              if parse_precision < 0:
                rhs=parsed_value.split(".")[1]
                parse_precision=len(rhs)

              #Calculate tolerance
              if parse_tolerance is not None:
                if parse_tolerance < 0:
                  lhs=parsed_value.replace("-","").split(".")[0]
                  rhs=parsed_value.replace("-","").split(".")[1]
                  if len(lhs) > 4:
                    parse_tolerance=5.0*10**(len(lhs)-12)
                  else:
                    parse_tolerance=5.0*10**(-(len(rhs)-2))

              parsed_values.append(
                create_JSON(
                  value_name=parse_name_line,
                  value_value=parsed_value,
                  value_header=parse_header,
                  value_type=parse_type,
                  value_precision=parse_precision,
                  value_tolerance=parse_tolerance))
            else:
              #Break out of the for loop over parse instances
              break

          if parsed_value is None:
            if parse_required:
              #Break out of the for loop over parse group
              break

        elif parse_multiple_indexes is not None:
          for parse_value_index in range(0,parse_multiple_indexes+1):
            parse_name_line=parse_name+" "+parse_multiple_label+" "+str(parse_value_index+1)
            file_handle.seek(0)
            parsed_value=parser(
              regex_string,file_handle,
              parse_memory_map=parse_memory_map_temp,
              parse_value_index=parse_value_index,
              parse_line_offset=parse_line_offset,
              parse_instance=parse_instance,
              parse_integer=parse_integer,
              parse_case_sensitive=parse_case_sensitive,
              parse_multiline=parse_multiline)

            if parsed_value is not None:

              #Calculate precision
              if parse_precision < 0:
                rhs=parsed_value.split(".")[1]
                parse_precision=len(rhs)

              #Calculate tolerance
              if parse_tolerance is not None:
                if parse_tolerance < 0:
                  lhs=parsed_value.replace("-","").split(".")[0]
                  rhs=parsed_value.replace("-","").split(".")[1]
                  if len(lhs) > 4:
                    parse_tolerance=5.0*10**(len(lhs)-12)
                  else:
                    parse_tolerance=5.0*10**(-(len(rhs)-2))

              parsed_values.append(
                create_JSON(
                  value_name=parse_name_line,
                  value_value=parsed_value,
                  value_header=parse_header,
                  value_type=parse_type,
                  value_precision=parse_precision,
                  value_tolerance=parse_tolerance))
            else:
              #Break out of the for loop over parse instances
              break

          if parsed_value is None:
            if parse_required:
              #Break out of the for loop over parse group
              break

      else:
        file_handle.seek(0)
        parsed_value=parser(
          regex_string,file_handle,
          parse_memory_map=parse_memory_map_temp,
          parse_value_index=parse_value_index,
          parse_line_offset=parse_line_offset,
          parse_instance=parse_instance,
          parse_integer=parse_integer,
          parse_case_sensitive=parse_case_sensitive,
          parse_multiline=parse_multiline)

        if parsed_value is not None:

          #Calculate precision
          if parse_precision < 0:
            rhs=parsed_value.split(".")[1]
            parse_precision=len(rhs)

          #Calculate tolerance
          if parse_tolerance is not None:

            if parse_tolerance < 0:

              #If we have an array
              if hasattr(parsed_value, '__len__') and (not isinstance(parsed_value, str)):
                lhs=parsed_value[0].replace("-","").split(".")[0]
                rhs=parsed_value[0].replace("-","").split(".")[1]
              #Else we have a scalar
              else:
                lhs=parsed_value.replace("-","").split(".")[0]
                rhs=parsed_value.replace("-","").split(".")[1]

              if len(lhs) > 4:
                parse_tolerance=5.0*10**(len(lhs)-12)
              else:
                parse_tolerance=5.0*10**(-(len(rhs)-2))

          parsed_values.append(
            create_JSON(
              value_name=parse_name,
              value_value=parsed_value,
              value_header=parse_header,
              value_type=parse_type,
              value_precision=parse_precision,
              value_tolerance=parse_tolerance))

        else:
          if parse_required:
            #Break out of the for loop over parse group
            break

  return parsed_values
