function compare_all_texts_in_test_folder( foldername)

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

% Delete all items except .txt files from list
i = 1;
while i <= length( henk)
  if endsWith( henk(i).name, '.txt')
    i = i+1;
  else
    henk( i) = [];
  end
end

% Compare all .txt files
all_files_match = true;
for i = 1: length( henk)

  disp(['Comparing text file ' num2str(i) '/' num2str( length( henk)) ': ' henk(i).name])

  filename_ref = [foldername_ref '/' henk(i).name];
  filename_mod = [foldername_mod '/' henk(i).name];

  if ~exist( filename_mod,'file')
    error(['file "' henk(i).name '" does not exist in results of test "' foldername '"'])
  end

  files_match = compare_text_file( filename_ref, filename_mod);

  if ~files_match
    all_files_match = false;
  end

end

if ~all_files_match
  error('Not all files are identical')
end

end