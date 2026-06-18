clc
clear all
close all

% Read test results
addpath(['/Users/Beren017/Documents/GitHub/UPSY-models/' ...
  'automated_testing/UPSY/component_test_PETSc_matrix_solving/results/'])
[KSPs, PCs, matrix_names, results_raw] = PETSc_solver_results;

% Process results
results_processed = process_results( KSPs, PCs, matrix_names, results_raw);

% Plot results
plot_results( KSPs, PCs, results_processed)

function results_processed = process_results( KSPs, PCs, matrix_names, results_raw)

nmats = size( matrix_names,1);
nKSPs = size( KSPs,1);
nPCs  = size( PCs ,1);

% Structure the data into nice arrays

n_Axb_its = zeros( nKSPs, nPCs, nmats);
rmse      = zeros( nKSPs, nPCs, nmats);
tcomp     = zeros( nKSPs, nPCs, nmats);

for i = 1: length( results_raw)

  mi   = find( strcmpi( matrix_names, results_raw(i).matrix_name));
  kspi = find( strcmpi( KSPs        , results_raw(i).KSP));
  pci  = find( strcmpi( PCs         , results_raw(i).PC));

  n_Axb_its( kspi, pci, mi) = results_raw(i).n_Axb_its;
  rmse(      kspi, pci, mi) = results_raw(i).rmse;
  tcomp(     kspi, pci, mi) = results_raw(i).tcomp;

end

[n_Axb_its_score_min_rel, n_Axb_its_score_max_rel] = ...
  calc_scores_min_max_rel( n_Axb_its);
[rmse_score_min_rel, rmse_score_max_rel] = ...
  calc_scores_min_max_rel( rmse);
[tcomp_score_min_rel, tcomp_score_max_rel] = ...
  calc_scores_min_max_rel( tcomp);

% Gather processed results for output
results_processed.n_Axb_its_score_min_rel = n_Axb_its_score_min_rel;
results_processed.n_Axb_its_score_max_rel = n_Axb_its_score_max_rel;
results_processed.rmse_score_min_rel      = rmse_score_min_rel;
results_processed.rmse_score_max_rel      = rmse_score_max_rel;
results_processed.tcomp_score_min_rel     = tcomp_score_min_rel;
results_processed.tcomp_score_max_rel     = tcomp_score_max_rel;

end

function [d_score_min_rel, d_score_max_rel] = calc_scores_min_max_rel( d)

nKSPs = size( d,1);
nPCs  = size( d,2);
nmats = size( d,3);

% For every matrix equation, rank all the KSP/PC combi's
% in terms of their number of iterations, and their rmse
% (lower = better)

d_score = zeros( nKSPs, nPCs, nmats);

for mi = 1: nmats
  d_score(:,:,mi) = calc_scores( d( :,:,mi));
end

% Reduce this to the best and worst scores for each KSP/PC
% combi across all matrices

d_score_min = zeros( nKSPs, nPCs);
d_score_max = zeros( nKSPs, nPCs);

for kspi = 1: nKSPs
  for pci = 1: nPCs
    vals = d_score( kspi, pci, :);
    d_score_min( kspi, pci) = min( vals(:));
    d_score_max( kspi, pci) = max( vals(:));
  end
end

% Scale scores to [0,1]
d_score_min_rel = zeros( nKSPs, nPCs);
d_score_max_rel = zeros( nKSPs, nPCs);

for kspi = 1: nKSPs
  for pci = 1: nPCs
    d_score_min_rel( kspi, pci) = ...
      (d_score_min( kspi, pci) - min( d_score_min (:))) / ...
      (max( d_score_min (:)) - min( d_score_min (:)));
    d_score_max_rel( kspi, pci) = ...
      (d_score_max( kspi, pci) - min( d_score_max (:))) / ...
      (max( d_score_max (:)) - min( d_score_max (:)));
  end
end

end

function scores = calc_scores( vals)
% Assign scores, including 'shared' places

vals_vec = vals(:);

[vals_sorted,ind] = sortrows( vals_vec);

scores_sorted = zeros( size( vals_vec));
score = 1;
scores_sorted( 1) = 1;
for i = 2: length( vals_vec)
  if vals_sorted( i) > vals_sorted( i-1)
    score = score+1;
  end
  scores_sorted( i) = score;
end

scores_vec = scores_sorted*0;
for i = 1: length( scores_sorted)
  ii = ind( i);
  scores_vec( ii) = scores_sorted( i);
end

scores = reshape( scores_vec, size( vals));

end

function plot_results( KSPs, PCs, results_processed)

nKSPs = size( KSPs,1);
nPCs  = size( PCs ,1);

wa = [400,200];
ha = 400;

margins_hor = [120,50,25];
margins_ver = [25,70];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);

ax1 = H.Ax{ 1,1};
ax2 = H.Ax{ 1,2};

set( ax2,'color','none')

pos = get( ax1,'position');
H.Colorbar = colorbar( ax1,'location','eastoutside');
set( ax1,'position',pos);

set( ax1,'xlim',[0,nPCs],'ylim',[0,nKSPs],'xgrid','off','ygrid','off',...
  'xtick',(0:nPCs -1)+0.5,'xticklabels',PCs,...
  'ytick',(0:nKSPs-1)+0.5,'yticklabels',KSPs);

ncols = 255;
cmap = crameri('vik',ncols);
cmap = [cmap(:,1), cmap(:,3), cmap(:,2)];
colormap( ax1, flipud( cmap));
set( H.Colorbar,'ytick',[0.1,0.9],'ticklabels',{'Bad','Good'});

markersize = 15;

for kspi = 1: nKSPs
  for pci = 1: nPCs

    % n_Axb_its_max (left)
    score = results_processed.n_Axb_its_score_max_rel( kspi, pci);
    x = pci  - 0.67;
    y = kspi - 0.33;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',ax1,'xdata',x,'ydata',y,'markersize',markersize,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

    % tcomp_max (right)
    score = results_processed.tcomp_score_max_rel( kspi, pci);
    x = pci  - 0.33;
    y = kspi - 0.33;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',ax1,'xdata',x,'ydata',y,'markersize',markersize,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

    % rmse_max (bottom)
    score = results_processed.rmse_score_max_rel( kspi, pci);
    x = pci  - 0.5;
    y = kspi - 0.67;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',ax1,'xdata',x,'ydata',y,'markersize',markersize+5,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

  end
end

%% Legend
pos = get( ax2,'position');
wa = pos( 3);
ha = pos( 4);
rxy = ha / wa;
set( ax2,'xlim',[-1,1],'ylim',[-rxy,rxy],...
  'xtick',[],'ytick',[]);
ax2.XAxis.Visible = 'off';
ax2.YAxis.Visible = 'off';

theta = linspace( 0, 2*pi, 100);
r0 = 0.7;
r1 = 0.35;
theta1 = 0.7 * pi;
theta2 = 0.3 * pi;
lw = 1.5;

% rmse_max (bottom)
x0 = 0;
y0 = -0.3;
text(ax2,x0-0.35,y0,'RMSE','fontsize',24)
xx = x0 + r0 * cos( theta);
yy = y0 + r0 * sin( theta);
patch('parent',ax2,'xdata',xx,'ydata',yy,'facecolor','none',...
  'edgecolor','k','linewidth',lw)

% n_Axb_its_max (left)
x1 = x0 + (r0+r1) * cos( theta1);
y1 = y0 + (r0+r1) * sin( theta1);
text(ax2,x1-0.25,y1,'# its.','fontsize',24)
xx = x1 + r1 * cos( theta);
yy = y1 + r1 * sin( theta);
patch('parent',ax2,'xdata',xx,'ydata',yy,'facecolor','none',...
  'edgecolor','k','linewidth',lw)

% tcomp_max (right)
r2 = r1;
x2 = x0 + (r0+r2) * cos( theta2);
y2 = y0 + (r0+r2) * sin( theta2);
text(ax2,x2-0.28,y2,'t_{comp}','fontsize',24)
xx = x2 + r2 * cos( theta);
yy = y2 + r2 * sin( theta);
patch('parent',ax2,'xdata',xx,'ydata',yy,'facecolor','none',...
  'edgecolor','k','linewidth',lw)

end


