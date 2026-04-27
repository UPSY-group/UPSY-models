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


henk = dir( foldername_ref);

% Delete all items except checksum netcdf files from list
i = 1;
while i <= length( henk)
  if endsWith( henk(i).name, '_checksum.nc')
    i = i+1;
  else
    henk( i) = [];
  end
end

% Compare all checksum netcdf files
all_files_match = true;
for i = 1: length( henk)

  disp(['Comparing netcdf file ' num2str(i) '/' num2str( length( henk)) ': ' henk(i).name])

  filename_ref = [foldername_ref '/' henk(i).name];
  filename_mod = [foldername_mod '/' henk(i).name];

  if ~exist( filename_mod,'file')
    error(['file "' henk(i).name '" does not exist in results of test "' foldername '"'])
  end

  files_match = compare_netcdf( filename_ref, filename_mod);

  if ~files_match
    all_files_match = false;
  end
  
end

if ~all_files_match
  error('Not all files are identical')
end

end