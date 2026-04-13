function reduce_netcdf_to_checksum( filename)
% Reduce data in a NetCDF file to only checksums

f_raw = ncinfo( filename);
f = f_raw;

[f, dim_checksum] = add_dim_checksum( f);

f.Variables = [];
f.Groups = [];
[f, original_names] = replace_variables_with_checksums( '', [], f_raw, f, dim_checksum);
filename_redux = create_checksum_file( filename, f);
checksum_all_variables( filename, filename_redux, f, original_names);

function [f,dim_checksum] = add_dim_checksum( f)
  dim_checksum.Name      = 'checksum';
  dim_checksum.Length    = 4;
  dim_checksum.Unlimited = false;
  f.Dimensions( end+1) = dim_checksum;
end

function [f, original_names] = replace_variables_with_checksums( ...
    parent_name, original_names, f_raw, f, dim_checksum)

  for vari = 1: length( f_raw.Variables)
  
    var = f_raw.Variables( vari);

    if ~isempty( parent_name)
      original_name = [parent_name '/' var.Name];
      var_name_new  = [parent_name '_' var.Name];
    else
      original_name = var.Name;
      var_name_new  = var.Name;
    end
  
    var.Name         = var_name_new;
    var.Dimensions   = dim_checksum;
    var.Size         = 4;
    var.Datatype     = 'double';
    var.FillValue    = 9e9;
    var.ChunkSize    = var.Size;
    var.Shuffle      = false;
    var.DeflateLevel = [];
  
    if isempty( f.Variables)
      f.Variables = var;
      original_names{1} = original_name;
    else
      f.Variables(    end+1) = var;
      original_names{ end+1} = original_name;
    end
  
  end

  % Then grouped variables - convert these into
  % regular ones, as Matlab has a bug in the group code...

  for grpi = 1: length( f_raw.Groups)
    grp = f_raw.Groups( grpi);
    if isempty( parent_name)
      group_name = grp.Name;
    else
      group_name = [parent_name '_' grp.Name];
    end
    [f, original_names] = replace_variables_with_checksums( ...
      group_name, original_names, grp, f, dim_checksum);
  end

end

function filename_redux = create_checksum_file( filename, f)

  filename_redux = [filename( 1:end-3) '_checksum.nc'];
  f.Filename = filename_redux;
  
  if exist( filename_redux,'file')
    delete( filename_redux)
  end
  
  ncwriteschema( filename_redux, f);

end

function checksum_all_variables( filename, filename_redux, f, original_names)
  
  for vari = 1: length( f.Variables)
  
    var = f.Variables( vari);
  
    % Read original variable
    d = ncread( filename, original_names{ vari});
    d = d(:);
  
    % Calculate checksums
    d_sum     = sum(      d);
    d_sum_abs = sum( abs( d));
    d_min     = min(      d);
    d_max     = max(      d);
  
    s = [d_sum, d_sum_abs, d_min, d_max];
  
    % Write to redux file
    ncwrite( filename_redux, var.Name, s);
  
  end

end

end