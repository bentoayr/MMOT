number_of_repetitions = 100;
numV = 16; 
points_per_hyper_edge = 3; 
num_hyp_edges = 500; 
num_points_per_class = 10;
num_classes = 7;
num_points = num_points_per_class*num_classes;
p = 0.05; 
fraction_viol = 0;  strength_of_viol = 0;  

all_errors_hyper = cell(number_of_repetitions,1);
all_errors_classic = cell(number_of_repetitions,1);

graph_vs_hypergraph = 3; 

list_of_hyper_edges_pairwise_MMOT = cell(number_of_repetitions,1);
list_of_hyper_edges_W_triang_area = cell(number_of_repetitions,1);
list_of_hyper_edges_W_barycenter = cell(number_of_repetitions,1);

thresholds_lev_range = [0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.3, 0.5, 0.98];


parfor rep_ix = 1:number_of_repetitions
    fprintf("Experiment %d out of %d\n",rep_ix,number_of_repetitions);

    rng('default'); 
    rng('shuffle');
    
    list_of_errors = cell(3,length(thresholds_lev_range));    
    
    [all_spectrums] =  generate_spectral_data(numV,num_points_per_class,num_classes,p);
    
    if (graph_vs_hypergraph == 1 || graph_vs_hypergraph == 3) 
        fprintf("Building hyperedges from MMOT that is an n-metric\n");
        hyper_edges_pairwise_MMOT = build_Hypergraph_edges(fraction_viol,strength_of_viol,all_spectrums,points_per_hyper_edge,num_hyp_edges,1);
        fprintf("Building hyperedges from MMOT that is not an n-metric \n");
        hyper_edges_W_triang_area = build_Hypergraph_edges(fraction_viol,strength_of_viol,all_spectrums,points_per_hyper_edge,num_hyp_edges,2);
        fprintf("Building hyperedges from MMOT that via barycenter\n");
        hyper_edges_W_barycenter = build_Hypergraph_edges(fraction_viol,strength_of_viol,all_spectrums,points_per_hyper_edge,num_hyp_edges,3);

        list_of_hyper_edges_pairwise_MMOT{rep_ix} = hyper_edges_pairwise_MMOT;
        list_of_hyper_edges_W_triang_area{rep_ix} = hyper_edges_W_triang_area;
        list_of_hyper_edges_W_barycenter{rep_ix} = hyper_edges_W_barycenter;
        
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
    
    if (graph_vs_hypergraph == 0 || graph_vs_hypergraph == 3)
       
        fprintf("Spectral clustering using WD\n");

        vec_data = [real(all_spectrums) ; imag(all_spectrums) ]'; 

        distances_data = 1000*ones(num_points);

        num_edges = num_hyp_edges*(3/2); 

        list_of_edges = nan(num_edges, 2); 
        for e =1:num_edges
            Result = 1;
            while Result
                g = randperm(num_points,2);
                [Result,~] = ismember(g,unique(list_of_edges(:,1:2),'rows'),'rows');
            end
            list_of_edges(e,1:2) = g;
        end
        
        for e = 1:num_edges
            i = list_of_edges(e,1);
            j = list_of_edges(e,2);
            distances_data(i,j) = compute_regular_WD_dist_uniform_points(vec_data(i,:),vec_data(j,:));
            distances_data(j,i) =  distances_data(i,j);
        end

        try
            spectral_not_reduced = spectralcluster((distances_data), num_classes, 'Distance', 'precomputed');
            label_err_spectral_not_reduced = compute_labeling_error(convert_index_to_matrix(spectral_not_reduced, num_classes),num_classes,num_points_per_class);
        catch
            label_err_spectral_not_reduced = 1 - 1/num_classes; % sometimes the method fails, in which case we output a random guess
        end

        all_errors_classic{rep_ix} = label_err_spectral_not_reduced;

    end
    
end

run('./generate_plots.m');