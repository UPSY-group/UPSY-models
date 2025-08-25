function H = plot_graph( graph)

% Axes and figure size
xw = graph.xmax - graph.xmin;
yw = graph.ymax - graph.ymin;
if xw >= yw
  wa = 800;
  ha = wa * yw / xw;
else
  ha = 800;
  wa = ha * xw / yw;
end
wf = 25 + wa + 100;
hf = 25 + ha + 50;

H.Fig = figure('position',[200,200,wf,hf],'color','w');
H.Ax  = axes('parent',H.Fig,'units','pixels','position',[25,25,wa,ha],'fontsize',24,...
  'xtick',[],'ytick',[],'xlim',[graph.xmin, graph.xmax],'ylim',[graph.ymin, graph.ymax]);

H.Graph = add_graph_to_plot( H.Ax, graph);

set( H.Ax,'units','normalized');

end