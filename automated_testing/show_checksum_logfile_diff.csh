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

set dir_ref = ${foldername}/${reference_dir}
set dir_mod = ${foldername}/results

# Check if expected directories exist
if (! -d $dir_ref) then
  echo "Error: Reference directory '$dir_ref' not found."
  exit 1
endif

if (! -d $dir_mod) then
  echo "Error: Results directory '$dir_mod' not found."
  exit 1
endif

# Gather all checksum logfiles in each directory.
# Supported pattern includes checksum_logfile.txt and checksum_logfile_*.txt.
set ref_files = (`find $dir_ref -maxdepth 1 -type f -name "checksum_logfile*.txt" -exec basename {} \; | sort`)
set mod_files = (`find $dir_mod -maxdepth 1 -type f -name "checksum_logfile*.txt" -exec basename {} \; | sort`)

if (($#ref_files == 0) && ($#mod_files == 0)) then
  echo "No checksum logfiles found in either '$dir_ref' or '$dir_mod'. Nothing to diff."
  exit 0
endif

# Require exact same set of checksum logfile names in both directories.
set sets_match = 1
if ($#ref_files != $#mod_files) then
  set sets_match = 0
else
  @ i = 1
  while ($i <= $#ref_files)
    if ("$ref_files[$i]" != "$mod_files[$i]") then
      set sets_match = 0
      break
    endif
    @ i ++
  end
endif

if (! $sets_match) then
  echo "Error: Reference and results folders do not contain the exact same set of checksum logfiles."
  echo "Reference files in '$dir_ref':"
  foreach f ($ref_files)
    echo "  $f"
  end
  echo "Results files in '$dir_mod':"
  foreach f ($mod_files)
    echo "  $f"
  end
  exit 1
endif

# Show diff for each matching checksum logfile.
foreach f ($ref_files)
  set filename_ref = ${dir_ref}/$f
  set filename_mod = ${dir_mod}/$f

  echo ""
  echo "=== Diff: $f ==="
  git diff --no-index $filename_ref $filename_mod
  set cmd_exit = $status

  if ($cmd_exit > 1) then
    echo "Error: git diff failed for '$f'."
    exit $cmd_exit
  endif
end

exit 0
endif
