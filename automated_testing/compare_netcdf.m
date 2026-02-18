function files_match = compare_netcdf( filename_ref, filename_mod)
% Compare two NetCDF files to see if they are identical,
% and if not, point out where the differences are

filename_ref_short = shorten_filename( filename_ref);
filename_mod_short = shorten_filename( filename_mod);
if ~strcmpi( filename_ref_short, filename_mod_short)
  disp(['Mismatching filenames: reference = "' ...
    filename_ref_short '", model = "' filename_mod_short '"'])
  files_match = false;
  return;
end

disp(['Comparing ' filename_ref_short])

info_matches = compare_ncinfo( filename_ref, filename_mod);
if info_matches
  files_match = compare_data( filename_ref, filename_mod);
else
  files_match = false;
end

  function filename_short = shorten_filename( filename)
    i = strfind( filename,'/');
    if isempty(i)
      filename_short = filename;
    else
      filename_short = filename( i(end)+1:end);
    end
  end

  function info_matches = compare_ncinfo( filename_ref, filename_mod)
    % Check if the header info of the two NetCDF files is identical

    f_ref = ncinfo( filename_ref);
    f_mod = ncinfo( filename_mod);

    info_matches = true;

    info_matches = info_matches && compare_ncinfo_attributes( ...
      shorten_filename( filename_ref), f_ref.Attributes, f_mod.Attributes);

    info_matches = info_matches && compare_ncinfo_dimensions( ...
      shorten_filename( filename_ref), f_ref.Dimensions, f_mod.Dimensions);

    info_matches = info_matches && compare_ncinfo_variables( ...
      shorten_filename( filename_ref), f_ref.Variables, f_mod.Variables);

    info_matches = info_matches && compare_ncinfo_groups( ...
      shorten_filename( filename_ref), f_ref.Groups, f_mod.Groups);

  end

  function info_matches = compare_ncinfo_dimensions( parent_name, dims_ref, dims_mod)
    % Check if the dimensions of the two NetCDF files are identical

    info_matches = true;

    if length( dims_ref) ~= length( dims_mod)
      disp('  Mismatching number of dimensions')
      info_matches = false;
      return
    end

    for dimi = 1: length( dims_ref)
      dim_ref = dims_ref( dimi);
      dim_mod = dims_mod( dimi);
      info_matches = info_matches && compare_dimensions( parent_name, dim_ref, dim_mod);
    end

  end

  function are_identical = compare_dimensions( parent_name, dim_ref, dim_mod)
    % Check if two NetCDF dimensions are identical

    are_identical = true;

    if ~strcmpi( dim_ref.Name, dim_mod.Name)
      are_identical = false;
      disp(['  Mismatching dimension name: reference = "' ...
        [parent_name '/' dim_ref.Name] '", model = "' [parent_name '/' dim_mod.Name] '"'])
    end

    if dim_ref.Length ~= dim_mod.Length
      are_identical = false;
      disp(['  Mismatching dimension "' [parent_name '/' dim_ref.Name] '" length: reference = ' ...
        num2str( dim_ref.Length) ', model = ' num2str( dim_mod.Length) ''])
    end

    if dim_ref.Unlimited ~= dim_mod.Unlimited
      are_identical = false;
      disp(['  Mismatching dimension "' [parent_name '/' dim_ref.Name] '" unlimitedness: reference = ' ...
        num2str( dim_ref.Unlimited) ', model = ' num2str( dim_mod.Unlimited) ''])
    end

  end

  function info_matches = compare_ncinfo_attributes( parent_name, atts_ref, atts_mod)
    % Check if the attributes of the two NetCDF files are identical

    info_matches = true;

    if length( atts_ref) ~= length( atts_mod)
      disp('  Mismatching number of attributes')
      info_matches = false;
      return
    end

    for atti = 1: length( atts_ref)
      att_ref = atts_ref( atti);
      att_mod = atts_mod( atti);
      info_matches = info_matches && compare_attributes( parent_name, att_ref, att_mod);
    end

  end

  function are_identical = compare_attributes( parent_name, att_ref, att_mod)
    % Check if two NetCDF attributes are identical

    are_identical = true;

    if ~strcmpi( att_ref.Name, att_mod.Name)
      are_identical = false;
      disp(['  Mismatching attribute name: reference = "' ...
        [parent_name '/' att_ref.Name] '", model = "' [parent_name '/' att_mod.Name] '"'])
      return
    end

    if ischar( att_ref.Value) && ischar( att_mod.Value)
      if ~strcmpi( att_ref.Value, att_mod.Value)
        are_identical = false;
        disp(['  Mismatching attribute "' [parent_name '/' att_ref.Name] '" value: reference = "' ...
          att_ref.Value '", model = "' att_mod.Value '"'])
      end
    else
      if att_ref.Value ~= att_mod.Value
        are_identical = false;
        disp(['  Mismatching attribute "' att_ref.Name '" value: reference = ' ...
          num2str( att_ref.Value) ', model = ' num2str( att_mod.Value)])
      end
    end

  end

  function info_matches = compare_ncinfo_variables( parent_name, vars_ref, vars_mod)
    % Check if the dimensions of the two NetCDF files are identical

    info_matches = true;

    if length( vars_ref) ~= length( vars_mod)
      disp('  Mismatching number of variables')
      info_matches = false;
      return
    end

    for vari = 1: length( vars_ref)
      var_ref = vars_ref( vari);
      var_mod = vars_mod( vari);
      info_matches = info_matches && compare_variables( parent_name, var_ref, var_mod);
    end

  end

  function are_identical = compare_variables( parent_name, var_ref, var_mod)
    % Check if two NetCDF variables are identical

    are_identical = true;

    % Name
    if ~strcmpi( var_ref.Name, var_mod.Name)
      are_identical = false;
      disp(['  Mismatching variable name: reference = "' ...
        [parent_name '/' var_ref.Name] '", model = "' [parent_name '/' var_mod.Name] '"'])
    end

    % Dimensions
    are_identical = are_identical && compare_ncinfo_dimensions( ...
      [parent_name '/' var_ref.Name], var_ref.Dimensions, var_mod.Dimensions);

    if length( var_ref.Size) ~= length( var_mod.Size)
      % Rank
      are_identical = false;
      disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" size: reference = [' ...
        num2str( var_ref.Size) '], model = [' num2str( var_mod.Size) ']'])
    else
      if ~all( var_ref.Size == var_mod.Size)
        % Size
        are_identical = false;
        disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" size: reference = [' ...
          num2str( var_ref.Size) '], model = [' num2str( var_mod.Size) ']'])
      end
    end

    % Datatype
    if ~strcmpi( var_ref.Datatype, var_mod.Datatype)
      are_identical = false;
      disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" datatype: reference = "' ...
        [parent_name '/' var_ref.Datatype] '", model = "' [parent_name '/' var_mod.Datatype] '"'])
    end

    % Attributes
    are_identical = are_identical && compare_ncinfo_attributes( ...
      [parent_name '/' var_ref.Name], var_ref.Attributes, var_mod.Attributes);

    % ChunkSize
    if var_ref.ChunkSize ~= var_mod.ChunkSize
      are_identical = false;
      disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" ChunkSize: reference = ' ...
        num2str( var_ref.ChunkSize) ', model = ' num2str( var_ref.ChunkSize)])
    end

    % FillValue
    if var_ref.FillValue ~= var_mod.FillValue
      are_identical = false;
      disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" FillValue: reference = ' ...
        num2str( var_ref.FillValue) ', model = ' num2str( var_ref.FillValue)])
    end

    % DeflateLevel
    if var_ref.DeflateLevel ~= var_mod.DeflateLevel
      are_identical = false;
      disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" DeflateLevel: reference = ' ...
        num2str( var_ref.DeflateLevel) ', model = ' num2str( var_ref.DeflateLevel)])
    end

    % Shuffle
    if var_ref.Shuffle ~= var_mod.Shuffle
      are_identical = false;
      disp(['  Mismatching variable "' [parent_name '/' var_ref.Name] '" Shuffle: reference = ' ...
        num2str( var_ref.Shuffle) ', model = ' num2str( var_ref.Shuffle)])
    end

  end

  function info_matches = compare_ncinfo_groups( parent_name, grps_ref, grps_mod)
    % Check if the groups of the two NetCDF files are identical

    info_matches = true;

    if isempty( grps_ref) && isempty( grps_mod)
      return
    end

    if length( grps_ref) ~= length( grps_mod)
      disp('  Mismatching number of groups')
      info_matches = false;
      return
    end

    for grpi = 1: length( grps_ref)
      grp_ref = grps_ref( grpi);
      grp_mod = grps_mod( grpi);
      info_matches = info_matches && compare_groups( parent_name, grp_ref, grp_mod);
    end

  end

  function are_identical = compare_groups( parent_name, grp_ref, grp_mod)
    % Check if two NetCDF groups are identical

    are_identical = true;

    % Name
    if ~strcmpi( grp_ref.Name, grp_mod.Name)
      are_identical = false;
      disp(['  Mismatching group name: reference = "' ...
        [parent_name '/' grp_ref.Name] '", model = "' [parent_name '/' grp_mod.Name] '"'])
    end

    % Dimensions
    are_identical = are_identical && compare_ncinfo_dimensions( ...
      [parent_name '/' grp_ref.Name], grp_ref.Dimensions, grp_mod.Dimensions);

    % Variables
    are_identical = are_identical && compare_ncinfo_variables( ...
      [parent_name '/' grp_ref.Name], grp_ref.Variables, grp_mod.Variables);

    % Attributes
    are_identical = are_identical && compare_ncinfo_attributes( ...
      [parent_name '/' grp_ref.Name], grp_ref.Attributes, grp_mod.Attributes);

    % Groups
    are_identical = are_identical && compare_ncinfo_groups( ...
      [parent_name '/' grp_ref.Name], grp_ref.Groups, grp_mod.Groups);

  end

  function data_matches = compare_data( filename_ref, filename_mod)
    % Check if the data of the two NetCDF files are identical

    f_ref = ncinfo( filename_ref);

    data_matches = true;

    for vari = 1: length( f_ref.Variables)
      var_name = f_ref.Variables( vari).Name;
      data_matches = data_matches && compare_data_variable( filename_ref, filename_mod, var_name);
    end

    for grpi = 1: length( f_ref.Groups)
      grp = f_ref.Groups( grpi);
      for vari = 1: length( grp.Variables)
        var_name = [grp.Name '/' grp.Variables( vari).Name];
        data_matches = data_matches && compare_data_variable( filename_ref, filename_mod, var_name);
      end
    end

  end

  function data_matches = compare_data_variable( filename_ref, filename_mod, var_name)
    % Check if the data of a variable the two NetCDF files are identical

    d_ref = ncread( filename_ref, var_name);
    d_mod = ncread( filename_mod, var_name);

    d_ref = d_ref(:);
    d_mod = d_mod(:);

    data_matches = all( d_ref == d_mod);

    if ~data_matches
      disp(['  Mismatching data in ' filename_ref '/' var_name])
    end

  end

end