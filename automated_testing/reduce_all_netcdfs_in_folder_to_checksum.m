function reduce_all_netcdfs_in_folder_to_checksum( foldername)

henk = dir( foldername);

for i = 1: length( henk)
  filename = [foldername '/' henk(i).name];
  disp(['Reducing file ' num2str(i) '/' num2str(length(henk)) ': ' filename])
  if endsWith( filename,'.nc')
    reduce_netcdf_to_checksum( filename);
    delete( filename)
  end
end

end