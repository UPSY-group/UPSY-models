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

for i = 1: length( results_raw)

  mi   = find( strcmpi( matrix_names, results_raw(i).matrix_name));
  kspi = find( strcmpi( KSPs        , results_raw(i).KSP));
  pci  = find( strcmpi( PCs         , results_raw(i).PC));

  n_Axb_its( kspi, pci, mi) = results_raw(i).n_Axb_its;
  rmse(      kspi, pci, mi) = results_raw(i).rmse;

end

% For every matrix equation, rank all the KSP/PC combi's
% in terms of their number of iterations, and their rmse
% (lower = better)

n_Axb_its_score = zeros( nKSPs, nPCs, nmats);
rmse_score      = zeros( nKSPs, nPCs, nmats);

for mi = 1: nmats

  n_Axb_its_score(:,:,mi) = calc_scores( n_Axb_its( :,:,mi));
  rmse_score(     :,:,mi) = calc_scores( rmse(      :,:,mi));

end

% Reduce this to the best and worst scores for each KSP/PC
% combi across all matrices

n_Axb_its_score_min = zeros( nKSPs, nPCs);
n_Axb_its_score_max = zeros( nKSPs, nPCs);
rmse_score_min      = zeros( nKSPs, nPCs);
rmse_score_max      = zeros( nKSPs, nPCs);

for kspi = 1: nKSPs
  for pci = 1: nPCs
    % n_Axb_its
    vals = n_Axb_its_score( kspi, pci, :);
    n_Axb_its_score_min( kspi, pci) = min( vals(:));
    n_Axb_its_score_max( kspi, pci) = max( vals(:));
    % rmse
    vals = rmse_score( kspi, pci, :);
    rmse_score_min( kspi, pci) = min( vals(:));
    rmse_score_max( kspi, pci) = max( vals(:));
  end
end

% Scale scores to [0,1]
n_Axb_its_score_min_rel = zeros( nKSPs, nPCs);
n_Axb_its_score_max_rel = zeros( nKSPs, nPCs);
rmse_score_min_rel      = zeros( nKSPs, nPCs);
rmse_score_max_rel      = zeros( nKSPs, nPCs);

for kspi = 1: nKSPs
  for pci = 1: nPCs
    % n_Axb_its
    n_Axb_its_score_min_rel( kspi, pci) = ...
      (n_Axb_its_score_min( kspi, pci) - min( n_Axb_its_score_min (:))) / ...
      (max( n_Axb_its_score_min (:)) - min( n_Axb_its_score_min (:)));
    n_Axb_its_score_max_rel( kspi, pci) = ...
      (n_Axb_its_score_max( kspi, pci) - min( n_Axb_its_score_max (:))) / ...
      (max( n_Axb_its_score_max (:)) - min( n_Axb_its_score_max (:)));
    % rmse
    rmse_score_min_rel( kspi, pci) = ...
      (rmse_score_min( kspi, pci) - min( rmse_score_min (:))) / ...
      (max( rmse_score_min (:)) - min( rmse_score_min (:)));
    rmse_score_max_rel( kspi, pci) = ...
      (rmse_score_max( kspi, pci) - min( rmse_score_max (:))) / ...
      (max( rmse_score_max (:)) - min( rmse_score_max (:)));
  end
end

% Gather processed results for output
results_processed.n_Axb_its_score_min_rel = n_Axb_its_score_min_rel;
results_processed.n_Axb_its_score_max_rel = n_Axb_its_score_max_rel;
results_processed.rmse_score_min_rel      = rmse_score_min_rel;
results_processed.rmse_score_max_rel      = rmse_score_max_rel;

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

wa = 300;
ha = 300;

margins_hor = [120,90];
margins_ver = [25,70];

H = setup_multipanel_figure( wa, ha, margins_hor, margins_ver);
H.Ax = H.Ax{ 1,1};
pos = get( H.Ax,'position');
H.Colorbar = colorbar( H.Ax,'location','eastoutside');
set( H.Ax,'position',pos);

set( H.Ax,'xlim',[0,nPCs],'ylim',[0,nKSPs],'xgrid','off','ygrid','off',...
  'xtick',(0:nPCs -1)+0.5,'xticklabels',PCs,...
  'ytick',(0:nKSPs-1)+0.5,'yticklabels',KSPs);

ncols = 255;
cmap = crameri('vik',ncols);
cmap = [cmap(:,1), cmap(:,3), cmap(:,2)];
colormap( H.Ax, flipud( cmap));
set( H.Colorbar,'ytick',[0.1,0.9],'ticklabels',{'Bad','Good'});

markersize = 15;

for kspi = 1: nKSPs
  for pci = 1: nPCs

    % n_Axb_its_min (top left)
    score = results_processed.n_Axb_its_score_min_rel( kspi, pci);
    x = pci  - 0.67;
    y = kspi - 0.33;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',H.Ax,'xdata',x,'ydata',y,'markersize',markersize,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

    % n_Axb_its_min (top right)
    score = results_processed.n_Axb_its_score_max_rel( kspi, pci);
    x = pci  - 0.33;
    y = kspi - 0.33;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',H.Ax,'xdata',x,'ydata',y,'markersize',markersize,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

    % rmse_min (bottom left)
    score = results_processed.rmse_score_min_rel( kspi, pci);
    x = pci  - 0.67;
    y = kspi - 0.67;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',H.Ax,'xdata',x,'ydata',y,'markersize',markersize,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

    % rmse_min (bottom right)
    score = results_processed.rmse_score_max_rel( kspi, pci);
    x = pci  - 0.33;
    y = kspi - 0.67;
    ci = 1 + floor( score * (ncols-1));
    col = cmap( ci,:);
    line('parent',H.Ax,'xdata',x,'ydata',y,'markersize',markersize,...
      'markerfacecolor',col,'markeredgecolor','k','marker','o')

  end
end

end


