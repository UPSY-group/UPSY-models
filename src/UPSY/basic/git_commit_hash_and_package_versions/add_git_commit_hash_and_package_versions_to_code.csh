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

# Determine PETSc version
set petsc_version = "UNKNOWN"
which pkg-config >& /dev/null
if ($status == 0) then
  set tmp = `pkg-config --modversion PETSc`
  if ($status == 0 && "$tmp" != "") then
    set petsc_version = "$tmp"
  else
    set tmp = `pkg-config --modversion petsc`
    if ($status == 0 && "$tmp" != "") set petsc_version = "$tmp"
  endif
endif

# Determine NetCDF version
set netcdf_version = "UNKNOWN"
which nc-config >& /dev/null
if ($status == 0) then
  set tmp = `nc-config --version | sed -E 's/^netCDF[[:space:]]+//'`
  if ($status == 0 && "$tmp" != "") set netcdf_version = "$tmp"
endif

# Determine OpenMPI version
set openmpi_version = "UNKNOWN"
which mpirun >& /dev/null
if ($status == 0) then
  set tmp = `mpirun --version | head -n 1 | sed -E 's/^mpirun \(Open MPI\)[[:space:]]*//'`
  if ($status == 0 && "$tmp" != "") set openmpi_version = "$tmp"
else
  which mpifort >& /dev/null
  if ($status == 0) then
    set tmp = `mpifort --showme:version | head -n 1 | sed -E 's/^Open MPI[[:space:]]*//'`
    if ($status == 0 && "$tmp" != "") set openmpi_version = "$tmp"
  endif
endif

printf "Compiling the code from git commit hash: \033[95m%s\033[0m\n" "$git_hash"
if ($has_changes == ".true.") then
  printf "  \033[33mWARNING: Repository has uncommitted changes\033[0m\n"
endif
printf "Using the following package versions:\n"
printf "  PETSc version  : \033[95m%s\033[0m\n" "$petsc_version"
printf "  NetCDF version : \033[95m%s\033[0m\n" "$netcdf_version"
printf "  OpenMPI version: \033[95m%s\033[0m\n" "$openmpi_version"
echo ""

# Escape values for robust sed replacement
set escaped_hash    = `echo "$git_hash" | sed 's/[&|\\]/\\&/g'`
set escaped_petsc   = `echo "$petsc_version" | sed 's/[&|\\]/\\&/g'`
set escaped_netcdf  = `echo "$netcdf_version" | sed 's/[&|\\]/\\&/g'`
set escaped_openmpi = `echo "$openmpi_version" | sed 's/[&|\\]/\\&/g'`

# Update all parameters, independent of their previous value
sed -E -i.bak "s#(character\(len=\*\),[[:space:]]*parameter :: git_commit_hash[[:space:]]*=[[:space:]]*)'[^']*'#\1'$escaped_hash'#" "$fortran_file"

# Update the has_uncommitted_changes parameter
sed -E -i.bak "s#(logical,[[:space:]]*parameter :: has_uncommitted_changes[[:space:]]*=[[:space:]]*)\.(true|false)\.#\1$has_changes#" "$fortran_file"
sed -E -i.bak "s#(character\(len=\*\),[[:space:]]*parameter :: petsc_version[[:space:]]*=[[:space:]]*)'[^']*'#\1'$escaped_petsc'#" "$fortran_file"
sed -E -i.bak "s#(character\(len=\*\),[[:space:]]*parameter :: netcdf_version[[:space:]]*=[[:space:]]*)'[^']*'#\1'$escaped_netcdf'#" "$fortran_file"
sed -E -i.bak "s#(character\(len=\*\),[[:space:]]*parameter :: openmpi_version[[:space:]]*=[[:space:]]*)'[^']*'#\1'$escaped_openmpi'#" "$fortran_file"

# Remove backup file
rm -f "$fortran_file.bak"

exit 0