clc
clear all
close all

mesh_basename = 'mesh_Ant_uniform_3.0000E+05_m_nit_Lloyd_2';

filename_animation = [mesh_basename '_focussing_animation'];

% MP4 video parameters
framerate      = 5;
movie_filename = 'filename_animation';

% Create .mp4 file
if exist( [movie_filename '.avi'],'file')
  delete( [movie_filename '.avi'])
end
VW = VideoWriter( movie_filename,'Uncompressed AVI');
VW.FrameRate = framerate;
open( VW);

H.Fig = figure('units','pixels','position',[100,100,500,500],'color','w');
H.Ax = axes('parent',H.Fig,'units','pixels','position',[10,10,480,480],...
  'xtick',[],'ytick',[],'xgrid','off','ygrid','off');
H.Ax.XAxis.Visible = 'off';
H.Ax.YAxis.Visible = 'off';

H.Patch = patch('parent',H.Ax,'vertices',[],'faces',[],'facecolor','none',...
  'edgecolor','k','linewidth',1);

filenames = list_focussed_mesh_filenames( mesh_basename);
for i = 1: length( filenames)
  filename = filenames{ i};
  mesh.V   = ncread( filename,'V');
  mesh.Tri = ncread( filename,'Tri');
  set( H.Patch,'vertices',mesh.V,'faces',mesh.Tri);
  drawnow('update')
  
  % Add frame to mp4
  writeVideo( VW, getframe( H.Fig));

end

% Close video object
close( VW);

function filenames = list_focussed_mesh_filenames( mesh_basename)

henk = dir( 'results');

filenames = {};
for i = 1: length( henk)
  if startsWith( henk(i).name, mesh_basename)
    filenames{end+1} = ['results/' henk(i).name];
  end
end

end