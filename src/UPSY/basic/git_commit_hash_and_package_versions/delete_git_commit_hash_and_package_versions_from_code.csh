#! /bin/csh -f

set fortran_file = "src/UPSY/basic/git_commit_hash_and_package_versions/git_commit_hash_and_package_versions.f90"

# Check if the file exists
if (! -f "$fortran_file") then
  echo "Error: File $fortran_file not found!"
  exit 1
endif

# Use sed to replace the actual commit hash back to INVALID
# This pattern matches any quoted string in the git_commit_hash parameter
sed -E -i.bak "s/(character\(len=\*\),[[:space:]]*parameter :: git_commit_hash[[:space:]]*=[[:space:]]*)'[^']*'/\1'INVALID'/" "$fortran_file"

# Reset has_uncommitted_changes back to .false.
sed -E -i.bak "s/(logical,[[:space:]]*parameter :: has_uncommitted_changes[[:space:]]*=[[:space:]]*)\.(true|false)\./\1.false./" "$fortran_file"

# Reset package versions back to INVALID.
sed -E -i.bak "s/(character\(len=\*\),[[:space:]]*parameter :: petsc_version[[:space:]]*=[[:space:]]*)'[^']*'/\1'INVALID'/" "$fortran_file"
sed -E -i.bak "s/(character\(len=\*\),[[:space:]]*parameter :: netcdf_version[[:space:]]*=[[:space:]]*)'[^']*'/\1'INVALID'/" "$fortran_file"
sed -E -i.bak "s/(character\(len=\*\),[[:space:]]*parameter :: openmpi_version[[:space:]]*=[[:space:]]*)'[^']*'/\1'INVALID'/" "$fortran_file"
sed -E -i.bak "s/(character\(len=\*\),[[:space:]]*parameter :: compiler_version[[:space:]]*=[[:space:]]*)'[^']*'/\1'INVALID'/" "$fortran_file"

# Remove backup file
rm -f "$fortran_file.bak"

exit 0