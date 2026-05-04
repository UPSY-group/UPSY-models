#!/bin/csh

if ($#argv == 1) then
  set foldername = "$argv[1]"
else
  echo "Usage: $0 [foldername]"
  exit 1
endif

set filename_ref = ${foldername}/reference/checksum_logfile.txt
set filename_mod = ${foldername}/results/checksum_logfile.txt

# Check if the file exists
if ((-e $filename_ref) && (-e $filename_mod)) then
  git diff --no-index $filename_ref $filename_mod
endif
