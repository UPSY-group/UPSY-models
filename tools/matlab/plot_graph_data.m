function H = plot_graph_data( graph, d)

if size(d,2) > 1
  d = mean( d,2);
end

H = plot_graph( graph);

pos = get( H.Ax,'position');
H.Cbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);

ncols = 256;
cmap = parula(ncols);
clim = [min(d), max(d)];
if (clim(1)==clim(2))
  clim(1) = clim(1)-1;
  clim(2) = clim(2)+1;
end
set( H.Ax,'clim',clim)

%% Regular nodes

for ci = 1: ncols
  linedata(ci).x = zeros( graph.n,1);
  linedata(ci).y = zeros( graph.n,1);
  linedata(ci).n = 0;
end

for ni = 1: graph.n
  if ~graph.is_ghost( ni)
    if isnan( d( ni)); continue; end
    ci = (d( ni) - clim(1)) / (clim(2) - clim(1));
    ci = max(1,min(ncols,1 + round(ci*(ncols-1))));

    linedata(ci).n = linedata(ci).n + 1;
    linedata(ci).x( linedata(ci).n) = graph.V( ni,1);
    linedata(ci).y( linedata(ci).n) = graph.V( ni,2);
  end
end

for ci = 1: ncols
  linedata(ci).x = linedata(ci).x( 1: linedata(ci).n);
  linedata(ci).y = linedata(ci).y( 1: linedata(ci).n);
end

for ci = 1: ncols
  line('parent',H.Ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor','k','markersize',12);
end

%% Ghost nodes

for ci = 1: ncols
  linedata(ci).x = zeros( graph.n,1);
  linedata(ci).y = zeros( graph.n,1);
  linedata(ci).n = 0;
end

for ni = 1: graph.n
  if graph.is_ghost( ni)
    ci = (d( ni) - clim(1)) / (clim(2) - clim(1));
    ci = max(1,min(ncols,1 + round(ci*(ncols-1))));

    linedata(ci).n = linedata(ci).n + 1;
    linedata(ci).x( linedata(ci).n) = graph.V( ni,1);
    linedata(ci).y( linedata(ci).n) = graph.V( ni,2);
  end
end

for ci = 1: ncols
  linedata(ci).x = linedata(ci).x( 1: linedata(ci).n);
  linedata(ci).y = linedata(ci).y( 1: linedata(ci).n);
end

for ci = 1: ncols
  line('parent',H.Ax,'xdata',linedata( ci).x,'ydata',linedata( ci).y,'linestyle','none',...
    'marker','o','markerfacecolor',cmap( ci,:),'markeredgecolor',cmap( ci,:),'markersize',12);
end

end