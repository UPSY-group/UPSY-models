function H = add_graph_to_plot( ax, graph, varargin)

% Regular nodes
nn = sum( graph.nC( ~graph.is_ghost)) * 3;
xx = NaN( nn,1);
yy = NaN( nn,1);

n = 0;
for ni = 1: graph.n
  if (~graph.is_ghost( ni))
    for ci = 1: graph.nC( ni)
      nj = graph.C( ni,ci);
      if (~graph.is_ghost( nj) && ni > nj)
        xx( n+1:n+3) = [graph.V( ni,1); graph.V( nj,1); NaN];
        yy( n+1:n+3) = [graph.V( ni,2); graph.V( nj,2); NaN];
        n = n+3;
      end
    end
  end
end

% Ghost nodes
ng = sum( graph.nC( graph.is_ghost)) * 3;
xxg = NaN( ng,1);
yyg = NaN( ng,1);

n = 0;
for ni = 1: graph.n
  if (graph.is_ghost( ni))
    for ci = 1: graph.nC( ni)
      nj = graph.C( ni,ci);
      if (graph.is_ghost( nj) && ni > nj) || (~graph.is_ghost( nj))
        xxg( n+1:n+3) = [graph.V( ni,1); graph.V( nj,1); NaN];
        yyg( n+1:n+3) = [graph.V( ni,2); graph.V( nj,2); NaN];
        n = n+3;
      end
    end
  end
end

%% Plot

linewidth       = 1;
color           = 'k';
markersize      = 5;

for vi = 1: 2: length( varargin)
  switch varargin{vi}
    case 'linewidth'
      linewidth = varargin{vi+1};
    case 'color'
      color = varargin{vi+1};
    case 'markersize'
      markersize = varargin{vi+1};
  end
end

H(1) = line('parent',ax,'xdata',xx,'ydata',yy,...
  'linewidth',linewidth,'linestyle','-','color',color,...
  'marker','o','markersize',markersize,...
  'markerfacecolor',color,'markeredgecolor',color);
H(2) = line('parent',ax,'xdata',xxg,'ydata',yyg,...
  'linewidth',linewidth,'linestyle','--','color',color,...
  'marker','none');
H(3) = line('parent',ax,'xdata',graph.V( graph.is_ghost,1),'ydata',graph.V( graph.is_ghost,2),...
  'linestyle','none',...
  'marker','o','markersize',markersize+3,...
  'markerfacecolor','none','markeredgecolor',color);

end