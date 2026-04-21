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
temp = textscan( fid,'%s','delimiter','\n', 'whitespace', '');
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


if     startsWith( strtrim( single_line),'fixed_output_dir_config')
  single_line = fixed_output_dir_config( opts);

elseif startsWith( strtrim( single_line),'xmin_ANT_config')
  single_line = xmin_ANT_config( opts);

elseif startsWith( strtrim( single_line),'xmax_ANT_config')
  single_line = xmax_ANT_config( opts);

elseif startsWith( strtrim( single_line),'ymin_ANT_config')
  single_line = ymin_ANT_config( opts);

elseif startsWith( strtrim( single_line),'ymax_ANT_config')
  single_line = ymax_ANT_config( opts);

elseif startsWith( strtrim( single_line),'dx_refgeo_init_idealised_config')
  single_line = dx_refgeo_init_idealised_config( opts);

elseif startsWith( strtrim( single_line),'dx_refgeo_PD_idealised_config')
  single_line = dx_refgeo_PD_idealised_config( opts);

elseif startsWith( strtrim( single_line),'dx_refgeo_GIAeq_idealised_config')
  single_line = dx_refgeo_GIAeq_idealised_config( opts);

elseif startsWith( strtrim( single_line),'refgeo_idealised_ISMIP_HOM_L_config')
  single_line = refgeo_idealised_ISMIP_HOM_L_config( opts);

elseif startsWith( strtrim( single_line),'maximum_resolution_uniform_config')
  single_line = maximum_resolution_uniform_config( opts);

elseif startsWith( strtrim( single_line),'dx_square_grid_smooth_ANT_config')
  single_line = dx_square_grid_smooth_ANT_config( opts);

elseif startsWith( strtrim( single_line),'choice_stress_balance_approximation_config')
  single_line = choice_stress_balance_approximation_config( opts);

elseif startsWith( strtrim( single_line),'dx_output_grid_ANT_config')
  single_line = dx_output_grid_ANT_config( opts);
  
elseif startsWith( strtrim( single_line),'transects_ANT_config')
  single_line = transects_ANT_config( opts);


end

end



function single_line = fixed_output_dir_config( opts)
  single_line = ['fixed_output_dir_config = "automated_testing/UPSY/integrated_test_ISMIP_HOM/' ...
    'results_ISMIP_HOM_' opts.experiment '_' opts.length_scale '_' opts.Stokes_approx '"'];
end

function single_line = xmin_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'xmin_ANT_config = -160000.0';
    case '80'
      single_line = 'xmin_ANT_config = -80000.0';
    case '40'
      single_line = 'xmin_ANT_config = -40000.0';
    case '20'
      single_line = 'xmin_ANT_config = -20000.0';
    case '10'
      single_line = 'xmin_ANT_config = -10000.0';
    case '5'
      single_line = 'xmin_ANT_config = -5000.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = xmax_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'xmax_ANT_config = 160000.0';
    case '80'
      single_line = 'xmax_ANT_config = 80000.0';
    case '40'
      single_line = 'xmax_ANT_config = 40000.0';
    case '20'
      single_line = 'xmax_ANT_config = 20000.0';
    case '10'
      single_line = 'xmax_ANT_config = 10000.0';
    case '5'
      single_line = 'xmax_ANT_config = 5000.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = ymin_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'ymin_ANT_config = -160000.0';
    case '80'
      single_line = 'ymin_ANT_config = -80000.0';
    case '40'
      single_line = 'ymin_ANT_config = -40000.0';
    case '20'
      single_line = 'ymin_ANT_config = -20000.0';
    case '10'
      single_line = 'ymin_ANT_config = -10000.0';
    case '5'
      single_line = 'ymin_ANT_config = -5000.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = ymax_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'ymax_ANT_config = 160000.0';
    case '80'
      single_line = 'ymax_ANT_config = 80000.0';
    case '40'
      single_line = 'ymax_ANT_config = 40000.0';
    case '20'
      single_line = 'ymax_ANT_config = 20000.0';
    case '10'
      single_line = 'ymax_ANT_config = 10000.0';
    case '5'
      single_line = 'ymax_ANT_config = 5000.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = dx_refgeo_init_idealised_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'dx_refgeo_init_idealised_config = 8000.0';
    case '80'
      single_line = 'dx_refgeo_init_idealised_config = 4000.0';
    case '40'
      single_line = 'dx_refgeo_init_idealised_config = 2000.0';
    case '20'
      single_line = 'dx_refgeo_init_idealised_config = 1000.0';
    case '10'
      single_line = 'dx_refgeo_init_idealised_config = 500.0';
    case '5'
      single_line = 'dx_refgeo_init_idealised_config = 250.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = dx_refgeo_PD_idealised_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'dx_refgeo_PD_idealised_config = 8000.0';
    case '80'
      single_line = 'dx_refgeo_PD_idealised_config = 4000.0';
    case '40'
      single_line = 'dx_refgeo_PD_idealised_config = 2000.0';
    case '20'
      single_line = 'dx_refgeo_PD_idealised_config = 1000.0';
    case '10'
      single_line = 'dx_refgeo_PD_idealised_config = 500.0';
    case '5'
      single_line = 'dx_refgeo_PD_idealised_config = 250.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = dx_refgeo_GIAeq_idealised_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'dx_refgeo_GIAeq_idealised_config = 8000.0';
    case '80'
      single_line = 'dx_refgeo_GIAeq_idealised_config = 4000.0';
    case '40'
      single_line = 'dx_refgeo_GIAeq_idealised_config = 2000.0';
    case '20'
      single_line = 'dx_refgeo_GIAeq_idealised_config = 1000.0';
    case '10'
      single_line = 'dx_refgeo_GIAeq_idealised_config = 500.0';
    case '5'
      single_line = 'dx_refgeo_GIAeq_idealised_config = 250.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = refgeo_idealised_ISMIP_HOM_L_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'refgeo_idealised_ISMIP_HOM_L_config = 160000.0';
    case '80'
      single_line = 'refgeo_idealised_ISMIP_HOM_L_config = 80000.0';
    case '40'
      single_line = 'refgeo_idealised_ISMIP_HOM_L_config = 40000.0';
    case '20'
      single_line = 'refgeo_idealised_ISMIP_HOM_L_config = 20000.0';
    case '10'
      single_line = 'refgeo_idealised_ISMIP_HOM_L_config = 10000.0';
    case '5'
      single_line = 'refgeo_idealised_ISMIP_HOM_L_config = 5000.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = maximum_resolution_uniform_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'maximum_resolution_uniform_config = 8000.0';
    case '80'
      single_line = 'maximum_resolution_uniform_config = 4000.0';
    case '40'
      single_line = 'maximum_resolution_uniform_config = 2000.0';
    case '20'
      single_line = 'maximum_resolution_uniform_config = 1000.0';
    case '10'
      single_line = 'maximum_resolution_uniform_config = 500.0';
    case '5'
      single_line = 'maximum_resolution_uniform_config = 250.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = dx_square_grid_smooth_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'dx_square_grid_smooth_ANT_config = 8000.0';
    case '80'
      single_line = 'dx_square_grid_smooth_ANT_config = 4000.0';
    case '40'
      single_line = 'dx_square_grid_smooth_ANT_config = 2000.0';
    case '20'
      single_line = 'dx_square_grid_smooth_ANT_config = 1000.0';
    case '10'
      single_line = 'dx_square_grid_smooth_ANT_config = 500.0';
    case '5'
      single_line = 'dx_square_grid_smooth_ANT_config = 250.0';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = choice_stress_balance_approximation_config( opts)
  switch opts.Stokes_approx
    case 'SIASSA'
      single_line = 'choice_stress_balance_approximation_config = "SIA/SSA"';
    case 'DIVA'
      single_line = 'choice_stress_balance_approximation_config = "DIVA"';
    case 'BPA'
      single_line = 'choice_stress_balance_approximation_config = "BPA"';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = dx_output_grid_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'dx_output_grid_ANT_config = 1000.0';
    case '80'
      single_line = 'dx_output_grid_ANT_config = 500.0';
    case '40'
      single_line = 'dx_output_grid_ANT_config = 250.0';
    case '20'
      single_line = 'dx_output_grid_ANT_config = 125.0';
    case '10'
      single_line = 'dx_output_grid_ANT_config = 62.5';
    case '5'
      single_line = 'dx_output_grid_ANT_config = 31.25';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end

function single_line = transects_ANT_config( opts)
  switch opts.length_scale
    case '160'
      single_line = 'transects_ANT_config = "ISMIP-HOM,dx=1000.0"';
    case '80'
      single_line = 'transects_ANT_config = "ISMIP-HOM,dx=500.0"';
    case '40'
      single_line = 'transects_ANT_config = "ISMIP-HOM,dx=250.0"';
    case '20'
      single_line = 'transects_ANT_config = "ISMIP-HOM,dx=125.0"';
    case '10'
      single_line = 'transects_ANT_config = "ISMIP-HOM,dx=62.5"';
    case '5'
      single_line = 'transects_ANT_config = "ISMIP-HOM,dx=31.25"';
    otherwise
      error(['invalid opts.length_scale ' opts.length_scale])
  end
end




