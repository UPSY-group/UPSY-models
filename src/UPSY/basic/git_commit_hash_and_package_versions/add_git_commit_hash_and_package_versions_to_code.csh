#! /bin/csh -f

set fortran_file = "src/UPSY/basic/git_commit_hash_and_package_versions/git_commit_hash_and_package_versions.f90"

# Check if the file exists
if (! -f "$fortran_file") then
  echo "Error: File $fortran_file not found!"
  exit 1
endif

# Get the current git commit hash
set git_hash = `git rev-parse HEAD`
if ($status != 0) then
  echo "Error: Could not get git commit hash. Make sure you're in a git repository."
  exit 1
endif

# Check for uncommitted changes
git diff-index --quiet HEAD
if ($status == 0) then
  set has_changes = ".false."
else
  set has_changes = ".true."
endif

echo "Adding git commit hash: $git_hash"
if ($has_changes == ".true.") then
  echo "WARNING: Repository has uncommitted changes"
endif

# Use sed to replace INVALID with the actual commit hash
# Need to escape special characters in sed
set escaped_hash = `echo "$git_hash" | sed 's/[&/\]/\\&/g'`
sed -i.bak "s/character(len=\*), parameter :: git_commit_hash         = 'INVALID'/character(len=*), parameter :: git_commit_hash         = '$escaped_hash'/" "$fortran_file"

# Update the has_uncommitted_changes parameter
sed -i.bak "s/logical,          parameter :: has_uncommitted_changes = \.false\./logical,          parameter :: has_uncommitted_changes = $has_changes/" "$fortran_file"

# Remove backup file
rm -f "$fortran_file.bak"

exit 0