function analyse_component_tests_mesh_focussing( foldername_automated_testing, do_print_figures)
% Analyse the results of all the mesh focussing component tests

disp('    Analysing mesh focussing component tests...')
disp('')

foldername_results = [foldername_automated_testing '/component_tests/results/mesh_focussing'];
foldername_figures = [foldername_automated_testing '/component_tests/figures'];

% List all the files
henk = dir( foldername_results);

% List all the unfocussed meshes
unfocussed_mesh_names = {};
for i = 1: length( henk)
  if startsWith( henk(i).name,'mesh_') && endsWith( henk(i).name,'.nc')
    ii = strfind( henk(i).name,'_r');
    mesh_name = henk(i).name(1:ii(end)-5);
    if isempty( unfocussed_mesh_names)
      unfocussed_mesh_names{1} = mesh_name;
    else
      if ~strcmpi( unfocussed_mesh_names{end},mesh_name)
        unfocussed_mesh_names{end+1} = mesh_name;
      end
    end
  end
end

H = setup_multipanel_figure( 800, 800, [5,5], [5,5]);
H.Ax = H.Ax{ 1,1};
set( H.Ax,'xtick',[],'ytick',[],'xgrid','off','ygrid','off')
H.Patch = patch('faces',[],'vertices',[],'facecolor','none');

for i = 1: length( unfocussed_mesh_names)
  unfocussed_mesh_name = unfocussed_mesh_names{i};
  
  % MP4 video parameters
  framerate      = 5;
  movie_filename = [foldername_figures '/mesh_focussing_' unfocussed_mesh_name];
  
  % Create .mp4 file
  if exist( [movie_filename '.mp4'],'file')
    delete( [movie_filename '.mp4'])
  end
  VW = VideoWriter( movie_filename,'Uncompressed AVI');
  VW.FrameRate = framerate;
  open( VW);

  for ii = 1: length( henk)
    if startsWith( henk(ii).name, unfocussed_mesh_name)
      focussed_mesh_name = henk(ii).name;
      filename = [foldername_results '/' focussed_mesh_name];
      mesh = read_mesh_from_file( filename);
      set( H.Patch,'faces',mesh.Tri,'vertices',mesh.V)
      set( H.Ax,'xlim',[mesh.xmin,mesh.xmax],'ylim',[mesh.ymin,mesh.ymax])
      drawnow('update');
      writeVideo( VW, getframe( H.Fig));
    end
  end

  % Close video object
  close( VW);

end
  
end