function H = add_graph_to_plot( ax, graph, varargin)

nn = sum( graph.nC)*3;
xx = zeros( nn,1);
yy = zeros( nn,1);

n = 0;
for ni = 1: graph.n
  for ci = 1: graph.nC( ni)
    nj = graph.C( ni,ci);
    xx( n+1:n+3) = [graph.V( ni,1); graph.V( nj,1); NaN];
    yy( n+1:n+3) = [graph.V( ni,2); graph.V( nj,2); NaN];
    n = n+3;
  end
end

linewidth       = 1;
linestyle       = '-';
color           = 'k';
marker          = 'none';
markersize      = 8;
markerfacecolor = 'none';
markeredgecolor = 'k';

for vi = 1: 2: length( varargin)
  switch varargin{vi}
    case 'linewidth'
      linewidth = varargin{vi+1};
    case 'linestyle'
      linestyle = varargin{vi+1};
    case 'color'
      color = varargin{vi+1};
    case 'marker'
      marker = varargin{vi+1};
    case 'markersize'
      markersize = varargin{vi+1};
    case 'markerfacecolor'
      markerfacecolor = varargin{vi+1};
    case 'markeredgecolor'
      markeredgecolor = varargin{vi+1};
  end
end

H = line('parent',ax,'xdata',xx,'ydata',yy,...
  'linewidth',linewidth,'linestyle',linestyle,'color',color,...
  'marker',marker','markersize',markersize,...
  'markerfacecolor',markerfacecolor,'markeredgecolor',markeredgecolor);

end