function reduce_all_netcdfs_in_folder_to_checksum( foldername)

henk = dir( foldername);

% Delete all items except non-checksum netcdf files from list
i = 1;
while i <= length( henk)
  if endsWith( henk(i).name, '.nc') && ~endsWith( henk(i).name, '_checksum.nc') && ...
      ~strcmpi( henk(i).name, 'resource_tracking.nc')
    i = i+1;
  else
    henk( i) = [];
  end
end

% Reduce all non-checksum netcdf files
for i = 1: length( henk)

  disp(['Reducing file ' num2str(i) '/' num2str(length(henk)) ': ' filename])

  filename = [foldername '/' henk(i).name];
  reduce_netcdf_to_checksum( filename);
  delete( filename)
  
end

end