function compare_all_netcdfs_in_test_folder( foldername)

if ~exist( foldername,'dir')
  error(['couldnt find test "' foldername '"'])
end

foldername_ref = [foldername '/reference'];
if ~exist( foldername_ref,'dir')
  error(['couldnt find reference for test "' foldername '"'])
end
foldername_mod = [foldername '/results'];
if ~exist( foldername_mod,'dir')
  error(['couldnt find results for test "' foldername '"'])
end

all_files_match = true;

henk = dir( foldername_ref);
for i = 1: length( henk)
  filename_ref = [foldername_ref '/' henk(i).name];
  if endsWith( filename_ref,'_checksum.nc')
    filename_mod = [foldername_mod '/' henk(i).name];
    if ~exist( filename_mod,'file')
      error(['file "' henk(i).name '" does not exist in results of test "' foldername '"'])
    end
    files_match = compare_netcdf( filename_ref, filename_mod);
    if ~files_match
      all_files_match = false;
    end
  end
end

if ~all_files_match
  error('Not all files are identical')
end

end