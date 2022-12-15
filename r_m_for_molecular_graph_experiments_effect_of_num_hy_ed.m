% make sure you run run_me_for_molecular_graph_experiments.m before running
% this file


hist_means = [];
hist_std = [];
for num_hyp_edges = 50:50:600

% maybe override some of the parameters
number_of_repetitions = 100;
numV = 16; 
points_per_hyper_edge = 3; 
num_points_per_class = 10;
num_classes = 5;
num_points = num_points_per_class*num_classes;
fraction_viol = 0; strength_of_viol = 0;    

all_errors_hyper = cell(number_of_repetitions,1);
all_errors_classic = cell(number_of_repetitions,1);

graph_vs_hypergraph = 3; 

prune_flag = 2;
ignore_weights_flag = 0;
truncate_number = 0;


for rep_ix = 1:number_of_repetitions
    fprintf("Experiment %d out of %d\n",rep_ix,number_of_repetitions);

    rng('default'); 
    rng('shuffle');
    
    %thresholds_lev_range = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.98];
    thresholds_lev_range = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98];
    %thresholds_lev_range = [0.05:0.01:0.98];
    
    list_of_errors = cell(3,length(thresholds_lev_range));    
    
      
    
    if (graph_vs_hypergraph == 1 || graph_vs_hypergraph == 3) 
        fprintf("Read hyperedges for MMOT that is an n-metric\n");
        hyper_edges_pairwise_MMOT = list_of_hyper_edges_pairwise_MMOT{rep_ix}(1:num_hyp_edges,:);
        
        fprintf("Read hyperedges for MMOT that is not an n-metric \n");
        hyper_edges_W_triang_area = list_of_hyper_edges_W_triang_area{rep_ix}(1:num_hyp_edges,:);
       
        fprintf("Read hyperedges for MMOT that via barycenter\n");
        hyper_edges_W_barycenter = list_of_hyper_edges_W_barycenter{rep_ix}(1:num_hyp_edges,:);
        
        for type_of_W_distance = 1:3

             for th_ix = 1:length(thresholds_lev_range)
                    thresholds_lev = thresholds_lev_range(th_ix);

                    switch type_of_W_distance
                        case 1
                            [clusters_nhcut,clusters_ttm] = Hypergraph_clustering_th(num_points,hyper_edges_pairwise_MMOT,thresholds_lev,num_classes);
                        case 2
                            [clusters_nhcut,clusters_ttm] = Hypergraph_clustering_th(num_points,hyper_edges_W_triang_area,thresholds_lev,num_classes);
                        case 3
                            [clusters_nhcut,clusters_ttm] = Hypergraph_clustering_th(num_points,hyper_edges_W_barycenter,thresholds_lev,num_classes);
                    end

                    label_err_nhcut = compute_labeling_error(clusters_nhcut,num_classes,num_points_per_class);
                    label_err_ttm = compute_labeling_error(clusters_ttm,num_classes,num_points_per_class);

                    list_of_errors{type_of_W_distance,th_ix} = [label_err_nhcut,label_err_ttm];

            end
        end
        all_errors_hyper{rep_ix} = list_of_errors;
    end
    

    
end

try
    run('./generate_plots.m'); 
catch
end

    hist_means = [hist_means; [n_metric_NH_cut_mean, non_n_metric_NH_cut_mean, barycenter_NH_cut_mean, n_metric_TTM_mean, non_n_metric_TTM_mean, barycenter_TTM_mean]];
    hist_std = [hist_std; [n_metric_NH_cut_std, non_n_metric_NH_cut_std, barycenter_NH_cut_std, n_metric_TTM_std, non_n_metric_TTM_std, barycenter_TTM_std]];

end

%%
figure;
hold on;
for l = 1:size(hist_means,2)
    errorbar( [50:50:500, 550, 600],hist_means(:,l),hist_std(:,l)/sqrt(number_of_repetitions))
end
ylabel('Fraction of miss-classified graphs');
xlabel('# of hyperedges for clustering via hypergraph cuts');
title('Clustering molecules with metrics and non-metrics');

legend({['pairwise NH-cut'],['non metric NH-cut'],['barycenter  NH-cut'],['pairwise TTM'],['non metric TTM'],['barycenter  TTM']},'Location','northwest');
box on;
set(gca,'fontname','times');
set(gca,'fontsize',14);