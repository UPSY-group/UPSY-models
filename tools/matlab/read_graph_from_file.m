function graph = read_graph_from_file( filename)

graph.xmin           = ncread( filename,'xmin');
graph.xmax           = ncread( filename,'xmax');
graph.ymin           = ncread( filename,'ymin');
graph.ymax           = ncread( filename,'ymax');

graph.V              = ncread( filename,'V');
graph.nC             = ncread( filename,'nC');
graph.C              = ncread( filename,'C');

graph.ni2vi          = ncread( filename,'ni2vi');
graph.ni2ti          = ncread( filename,'ni2ti');
graph.ni2ei          = ncread( filename,'ni2ei');

graph.vi2ni          = ncread( filename,'vi2ni');
graph.ti2ni          = ncread( filename,'ti2ni');
graph.ei2ni          = ncread( filename,'ei2ni');

graph.is_border      = ncread( filename,'is_border') == 1;
graph.border_nhat    = ncread( filename,'border_nhat');

graph.n              = size( graph.V,1);

end