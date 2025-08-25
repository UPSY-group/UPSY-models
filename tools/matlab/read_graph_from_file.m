function graph = read_graph_from_file( filename)

graph.xmin           = ncread( filename,'xmin');
graph.xmax           = ncread( filename,'xmax');
graph.ymin           = ncread( filename,'ymin');
graph.ymax           = ncread( filename,'ymax');

graph.V              = ncread( filename,'V');
graph.nC             = ncread( filename,'nC');
graph.C              = ncread( filename,'C');

graph.is_ghost       = ncread( filename,'is_ghost') == 1;
graph.ghost_nhat     = ncread( filename,'ghost_nhat');

graph.n              = size( graph.V,1);

end