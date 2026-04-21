clc
clear all
close all

filename_config_template = 'AA_config_template.cfg';
c = read_config_template( filename_config_template);

experiments    = {'A','B','C','D'};
length_scales  = {'160','80','40','20','10','5'};
Stokes_approxs = {'SIASSA','DIVA','BPA'};

for ei = 1: length( experiments)
  for li = 1: length( length_scales)
    for si = 1: length( Stokes_approxs)

      opts.experiment    = experiments{    ei};
      opts.length_scale  = length_scales{  li};
      opts.Stokes_approx = Stokes_approxs{ si};

      cc = setup_config( c, opts);

      config_filename = ['config_ISMIP_HOM' ...
        '_' opts.experiment ...
        '_' opts.length_scale ...
        '_' opts.Stokes_approx ...
        '.cfg'];

      write_config_to_file( cc, config_filename)

    end
  end
end

function c = read_config_template( filename)

if ~exist( filename,'file')
  error(['Config template "' filename '" not found'])
end

fid = fopen( filename);
temp = textscan( fid,'%s','delimiter','\n');
c = temp{1};
fclose( fid);

end

function c = setup_config( c, opts)

for li = 1: size( c,1)
  c{li} = setup_config_line( c{li}, opts);
end

end

function write_config_to_file( cc, config_filename)

if exist( config_filename,'file')
  delete( config_filename)
end

fid = fopen( config_filename,'w');
for li = 1: size( cc,1)
  fprintf( fid,'%s\n', cc{li});
end
fclose( fid);

end

%% Individual config variables
function single_line = setup_config_line( single_line, opts)


if     startsWith( single_line,'start_time_of_run_config')
  single_line = start_time_of_run_config( opts);


end

end



function single_line = start_time_of_run_config( opts)
  single_line = 'start_time_of_run_config = 0.0';
end
