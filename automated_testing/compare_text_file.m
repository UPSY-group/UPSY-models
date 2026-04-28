function files_match = compare_text_file( filename_ref, filename_mod)
% Compare two text files to see if they are identical,
% and if not, point out where the first difference is

files_match = true;

filename_ref_short = shorten_filename( filename_ref);
filename_mod_short = shorten_filename( filename_mod);
if ~strcmpi( filename_ref_short, filename_mod_short)
  disp(['Mismatching filenames: reference = "' ...
    filename_ref_short '", model = "' filename_mod_short '"'])
  files_match = false;
  return;
end

% Read reference file
fid_ref = fopen( filename_ref, 'r');
temp = textscan( fid_ref, '%s', 'delimiter', '\n', 'whitespace', '');
text_ref = temp{1};
fclose( fid_ref);

% Read model file
fid_mod = fopen( filename_mod, 'r');
temp = textscan( fid_mod, '%s', 'delimiter', '\n', 'whitespace', '');
text_mod = temp{1};
fclose( fid_mod);

if length( text_ref) ~= length( text_mod)
  files_match = false;
  disp(['  Mismatching number of lines in ' filename_ref_short])
end

% Compare line-by-line
for i = 1: min( length( text_ref), length( text_mod))
  line_ref = text_ref{ i};
  line_mod = text_mod{ i};
  lines_match = compare_lines( line_ref, line_mod, filename_ref_short);
  files_match = files_match && lines_match;
  if ~lines_match
    % Stop at the first mismatching line
    break
  end
end

function filename_short = shorten_filename( filename)
  i = strfind( filename,'/');
  if isempty(i)
    filename_short = filename;
  else
    filename_short = filename( i(end)+1:end);
  end
end

function lines_match = compare_lines( line_ref, line_mod, filename)

  % Exceptions for the header lines showing the program info:
  %   show any differences, but do not fail the test
  if  startsWith( strtrim( line_ref), 'Git commit hash') || ...
      startsWith( strtrim( line_ref), 'has uncommitted changes') || ...
      startsWith( strtrim( line_ref), 'PETSc version') || ...
      startsWith( strtrim( line_ref), 'NetCDF version') || ...
      startsWith( strtrim( line_ref), 'OpenMPI version') || ...
      startsWith( strtrim( line_ref), 'Compiler') || ...
      startsWith( strtrim( line_ref), 'Compiler flags')

    lines_match = true;

    if ~strcmp( line_ref, line_mod)
      disp(['Mismatching program info in ' filename])
      disp(['  Ref: ' line_ref])
      disp(['  Mod: ' line_mod])
    end

  % Default line comparison
  else

    lines_match = strcmp( line_ref, line_mod);
  
    if ~lines_match
      disp(['Mismatching lines in ' filename])
      disp(['  Ref: ' line_ref])
      disp(['  Mod: ' line_mod])
    end

  end

end

end