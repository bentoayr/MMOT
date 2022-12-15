function [all_spectrums_2] =  generate_spectral_data(numV,num_points_per_class,num_classes,p)

    all_spectrums_2 = zeros(2*numV,num_points_per_class*num_classes);
    
    centroid_graph = cell(num_classes,1);
    for type_of_graph = 1:num_classes 
        [~,  numE, numEline, Adj_G,~, ~, ~, ~, ~, ~, ~, ~ ] = generate_graph_data(numV, type_of_graph);
        centroid_graph{type_of_graph} = {numV,  numE, numEline, Adj_G} ;
    end
    
    for type_of_graph = 1:num_classes
        adj_center = centroid_graph{type_of_graph}{4};
        for point_id = 1:num_points_per_class

            noise = rand(numV) > p;
            noise = triu(noise,1) + triu(noise,1)';
            adj = 0.5*(2*(adj_center - 0.5).*(2*((noise) - 0.5)) + 1);
            adj = adj - diag(diag(adj));
          
            [approx_non_back_track_matrix] = compute_non_backtracking_matrix(numV,adj);
            
            spec_approx_non_back_track_matrix = sort(eig(approx_non_back_track_matrix)); 
            all_spectrums_2(:,point_id + num_points_per_class*(type_of_graph-1)) = spec_approx_non_back_track_matrix;
        end
    end
end