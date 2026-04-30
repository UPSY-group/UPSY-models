function files_match = compare_program_info( foldername)
% Compare program info (git commit hash, package versions, compiler, etc.
% between reference and test)

infos = {...
  'git commit hash'
  'PETSc version'
  'NetCDF version'
  'OpenMPI version'
  'compiler version'
  'compiler flags'};

program_info_mod = read_program_info( [foldername '/results'  ], infos);
program_info_ref = read_program_info( [foldername '/reference'], infos);
show_both_program_infos( foldername, program_info_ref, program_info_mod);

  function program_info = read_program_info( foldername, fieldnames)

    filename = find_first_netcdf_file_in_dir( foldername);

    for fi = 1: length( fieldnames)

      attname = fieldnames{ fi};

      % Older references might not have all the program info; allow this
      try
        attval = ncreadatt( filename, '/', attname);
      catch
        attval = '';
      end

      fieldname = strrep( attname, ' ', '_');
      program_info.(fieldname) = attval;

    end

  end

  function filename = find_first_netcdf_file_in_dir( foldername)

    filelist = dir( foldername);
    
    foundone = false;
    for i = 1: length( filelist)
      filename_short = filelist( i).name;
      if endsWith( filename_short,'.nc')
        foundone = true;
        filename = [foldername '/' filename_short];
        break
      end
    end
    if ~foundone
      error(['Could not find a NetCDF file in ' foldername]);
    end
  end

  function show_both_program_infos( foldername, program_info_ref, program_info_mod)

    fprintf('\n')
    fprintf(['Program info for ' foldername '\n'])
    fprintf('\n')

    fieldnames = fields( program_info_ref);
    for fi = 1: length( fieldnames)
      fieldname = fieldnames{ fi};
      field_ref = program_info_ref.(fieldname);
      field_mod = program_info_mod.(fieldname);

      disp([' ' fieldname])
      disp(['  Reference: ' field_ref])
      disp(['  Results  : ' field_mod])
    end

  end

end