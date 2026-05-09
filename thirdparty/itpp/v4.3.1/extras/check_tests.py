#!/usr/bin/env python

# File:   check_tests.py
# Brief:  Run each test and compare its output with the reference output.
# Author: Bogdan Cristea
#
# Usage: ./check_tests.py -r path_to_ref_files -w path_to_test_binaries
#
# This script runs all test binaries found in path_to_test_binaries and 
# compares their output with the corresponding reference output found in 
# path_to_ref_files.
#
# -------------------------------------------------------------------------
#
# Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
#
# This file is part of IT++ - a C++ library of mathematical, signal
# processing, speech processing, and communications classes and functions.
#
# IT++ is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along
# with IT++.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------

import os
import sys
from subprocess import call
from optparse import OptionParser
import platform

parser = OptionParser()
parser.add_option("-r", "--ref-path", dest="ref_path",
                      help="path where reference files (*.ref) are stored")
parser.add_option("-w", "--work-path", dest="path",
                      help="path where test binaries are stored")

(options, args) = parser.parse_args()
if (1 == len(sys.argv)):
  parser.print_help()
  sys.exit(1)

ref_path = options.ref_path
path =  options.path

def compare_files(file_name_ref, file_name):
  fd_ref = open(file_name_ref)
  fd = open(file_name)

  identical = 0
  for line_ref in fd_ref:
    line = fd.readline()
    if (line_ref != line):
      identical = 1
      break
  #check if the second file is bigger than the first one
  if 0 == identical:
    for line in fd:
      line_ref = fd_ref.readline()
      if (line != line_ref):
        identical = 1
        break

  fd_ref.close()
  fd.close()
  return identical

dirList=os.listdir(path)
i = 0
passed_tests = 0
failed_tests = 0
failed_test_name = []

fd_null = open(os.devnull, "w")

for fname in dirList:
  full_name = path+"/"+fname
  if os.access(full_name, os.X_OK) and (not os.path.isdir(full_name)):
    name_ext = os.path.splitext(fname)
    if (1 < len(name_ext)) and ("Windows" == platform.system()) and not(".exe" == name_ext[-1]):
      continue #skip non exe (Windows only)
    sys.stdout.write("# " + str(i) + " Running " + fname)
    sys.stdout.flush()
    i = i+1
    fd = open(full_name+".tmp", "w")
    call(full_name, stdout=fd, stderr=fd_null)
    fd.close()
#    out = call(["diff", full_name+".tmp", ref_path+"/"+fname+".ref"], stdout=fd_null)
    out = compare_files(full_name+".tmp", ref_path+"/"+name_ext[0]+".ref")
    if 0 == out:
      print(" PASS")
      passed_tests = passed_tests+1
    else:
      print(" FAILED")
      failed_tests = failed_tests+1
      failed_test_name.append(name_ext[0])
fd_null.close()
print("*"*60)
print("From " + str(i) + " tests " + str(passed_tests) + " have passed and " + str(failed_tests) + " have failed")
if 0 != failed_tests:
  print("Failed tests:")
  for i in range(0, len(failed_test_name)):
    print(" "+failed_test_name[i])
print("*"*60)
