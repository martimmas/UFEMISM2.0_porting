function graph = read_graph_from_file( filename)

graph.V              = ncread( filename,'V');
graph.nC             = ncread( filename,'nC');
graph.C              = ncread( filename,'C');

graph.is_ghost       = ncread( filename,'is_ghost') == 1;
graph.ghost_nhat     = ncread( filename,'ghost_nhat');

graph.n              = size( graph.V,1);
graph.ng             = sum( graph.is_ghost);
graph.nn             = graph.n - graph.ng;

end