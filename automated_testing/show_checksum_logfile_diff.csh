#!/bin/csh

if (($#argv == 1) || ($#argv == 2)) then
  set foldername = "$argv[1]"
  if ($#argv == 2) then
    set reference_dir = "$argv[2]"
  else
    set reference_dir = "reference"
  endif
else
  echo "Usage: $0 [foldername] [reference_dir]"
  exit 1
endif

set filename_ref = ${foldername}/${reference_dir}/checksum_logfile.txt
set filename_mod = ${foldername}/results/checksum_logfile.txt

# Check if the file exists
if ((-e $filename_ref) && (-e $filename_mod)) then
  git diff --no-index $filename_ref $filename_mod
endif
