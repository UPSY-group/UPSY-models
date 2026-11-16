clc
clear all
close all

fields = {
  'so'
  'thetao'
  };

years = 2015:2024;

foldername_base_src = ['/Users/Beren017/Documents/' ...
  'Downloads_from_Ghub'];

foldername_base_dst = 'AIS/CESM2-WACCM/ssp585/ocean';

for fi = 1: length( fields)
  
  field = fields{fi};

  for yi = 1: length( years)

    year = years( yi);
    year_str = num2str( year);

    disp(['Upscaling ' field ' for year ' year_str '...'])

    filename_src = [foldername_base_src '/' ...
      field '_AIS_CESM2-WACCM_ssp585_ocean_v3_2015-2024.nc'];
    filename_dst = [foldername_base_dst '/' field '/v3/' ...
      field '_AIS_CESM2-WACCM_ssp585_ocean_v3_2015-2024.nc'];

    x_src = double( ncread( filename_src,'x'));
    y_src = double( ncread( filename_src,'y'));
    z_src = double( ncread( filename_src,'z'));
    [xg_src, yg_src, zg_src] = meshgrid( x_src, y_src, z_src);
    nx_src = length( x_src);
    ny_src = length( y_src);
    nz_src = length( z_src);

    dx_dst = 40e3;
    x_dst = -3040e3: dx_dst: 3040e3;
    y_dst = -3040e3: dx_dst: 3040e3;
    z_dst = z_src;
    [xg_dst, yg_dst, zg_dst] = meshgrid( x_dst, y_dst, z_dst);
    nx_dst = length( x_dst);
    ny_dst = length( y_dst);
    nz_dst = length( z_dst);

    time = ncread( filename_src,'time');
    time_double = double( time);
    nt = length( time);

    %% Adapt f to new grid size

    f = ncinfo( filename_src);

    for di = 1: length( f.Dimensions)
      dim = f.Dimensions( di);
      if strcmpi( dim.Name,'x')
        dim.Length = nx_dst;
      elseif strcmpi( dim.Name,'y')
        dim.Length = ny_dst;
      elseif strcmpi( dim.Name,'z')
        dim.Length = nz_dst;
      end
      f.Dimensions(di) = dim;
    end

    for vi = 1: length( f.Variables)

      var = f.Variables(vi);

      for di = 1: length( var.Dimensions)
        dim = var.Dimensions( di);
        if strcmpi( dim.Name,'x')
          dim.Length = nx_dst;
        elseif strcmpi( dim.Name,'y')
          dim.Length = ny_dst;
        elseif strcmpi( dim.Name,'z')
          dim.Length = nz_dst;
        end
        var.Dimensions(di) = dim;
      end

      for si = 1: length( var.Size)
        if var.Size(si) == nx_src
          var.Size(si) = nx_dst;
        elseif var.Size(si) == ny_src
          var.Size(si) = ny_dst;
        elseif var.Size(si) == nz_src
          var.Size(si) = nz_dst;
        end
      end

      var.ChunkSize    = [];
      var.DeflateLevel = 1;
      var.Shuffle      = true;

      f.Variables(vi) = var;

    end

    %% Create new file
    if exist( filename_dst,'file')
      delete( filename_dst)
    end
    ncwriteschema( filename_dst, f);

    %% Read/upscale/write data
    ncwrite( filename_dst,'x'   ,x_dst);
    ncwrite( filename_dst,'y'   ,y_dst);
    ncwrite( filename_dst,'z'   ,z_dst);
    ncwrite( filename_dst,'time',time);

    d_src = ncread( filename_src, field);

    d_dst = zeros( nx_dst, ny_dst, nz_dst, length (years));
    for ti = 1: length( time)
      d_dst( :,:,:,ti) = interp3( xg_src, yg_src, zg_src, d_src( :,:,:,ti), xg_dst, yg_dst, zg_dst);
    end

    ncwrite( filename_dst, field, d_dst);
    
  end

end