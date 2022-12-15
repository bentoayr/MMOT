function [hyper_edges] = build_Hypergraph_edges(fraction_viol,strength_of_viol,all_data_points,num_points_in_hyper_edge,num_hyp_edges,distance_func_type)

m = size(all_data_points,1);
num_points = size(all_data_points,2);

hyper_edges = nan(num_hyp_edges,num_points_in_hyper_edge+1);

for t=1: num_hyp_edges - 4*fraction_viol*num_hyp_edges
    
    fprintf(" %d%%| ", round(100*t/num_hyp_edges));
    
    Result = 1;
    while Result
        g = randperm(num_points,num_points_in_hyper_edge);
        [Result,~] = ismember(g,unique(hyper_edges(:,1:num_points_in_hyper_edge),'rows'),'rows');
    end
    hyper_edges(t,1:num_points_in_hyper_edge) = g;
    
    subpoints = nan(2, num_points_in_hyper_edge, m);
    for ii = 1:num_points_in_hyper_edge
        subpoints(1,ii,:) = real(all_data_points(:,g(ii)));
        subpoints(2,ii,:) = imag(all_data_points(:,g(ii)));
    end
    
    dist_BaryC = nan( m , m , m );
    dist_3 = nan( m , m , m );
    dist_2 = nan( m , m, m);
    dist_pairs = cell(3,1);
    dist_pairs{1} = nan(m,m); dist_pairs{2} = nan(m,m); dist_pairs{3} = nan(m,m);
    for i = 1:m
        for j = 1:m
            for k = 1:m
                dist_3(i,j,k) = compute_d_via_triangle_area([ subpoints( : , 1 , i ), subpoints( : , 2 , j ), subpoints( : , 3 , k ) ], eps);
                dist_2(i,j,k) = compute_d_via_pairwise_sums([ subpoints( : , 1 , i ), subpoints( : , 2 , j ), subpoints( : , 3 , k ) ]);
                dist_BaryC(i,j,k) = compute_d_via_barycenter([ subpoints( : , 1 , i ), subpoints( : , 2 , j ), subpoints( : , 3 , k ) ]);
                dist_pairs{1}(i,j) = norm(subpoints( : , 1 , i ) -  subpoints( : , 2 , j ));
                dist_pairs{2}(i,k) = norm(subpoints( : , 1 , i ) -  subpoints( : , 3 , k ));
                dist_pairs{3}(j,k) = norm(subpoints( : , 2 , j ) -  subpoints( : , 3 , k ));
            end
        end
    end
    
    try
        switch (distance_func_type)
            case 1
                probs = ones(3,m)/m;
                [W, ~] = compute_W_n_dist_for_generic_dist_func(probs,dist_2);
            case 2
                probs = ones(3,m)/m;
                [W, ~] = compute_W_n_dist_for_generic_dist_func(probs,dist_3);
            case 3
                probs = ones(3,m)/m;
                [W, ~] = compute_W_n_dist_for_generic_dist_func(probs,dist_BaryC);
        end
    catch
        W = 10000;  % very rare, ignore hyperedge
    end
    
    if (isnan(W)) % very rare, ignore hyperedge
        W = 10000;
    end
    hyper_edges(t,4) = W;
end

for t= 1 + num_hyp_edges - 4*fraction_viol*num_hyp_edges :4:num_hyp_edges
  fprintf(" %d%%| ", round(100*t/num_hyp_edges));

    Result = 1;
    while Result
        g = randperm(num_points,num_points_in_hyper_edge+1);
        [Result,~] = ismember(g,unique(hyper_edges(:,1:num_points_in_hyper_edge+1),'rows'),'rows');
    end
    
    sub_id_list = [[1,2,3];[1,2,4];[1,3,4];[2,3,4]]';
    for sub_t = 1:4
        sub_id = sub_id_list(:,sub_t);
        g_sub = g(sub_id);
        hyper_edges(t+sub_t-1,1:num_points_in_hyper_edge) = g_sub;
        
        subpoints = nan(2, num_points_in_hyper_edge, m);
        for ii = 1:num_points_in_hyper_edge
            subpoints(1,ii,:) = real(all_data_points(:,g_sub(ii)));
            subpoints(2,ii,:) = imag(all_data_points(:,g_sub(ii)));
        end

        dist_BaryC = nan( m , m , m );
        dist_3 = nan( m , m , m );
        dist_2 = nan( m , m, m);
        dist_pairs = cell(3,1);
        dist_pairs{1} = nan(m,m); dist_pairs{2} = nan(m,m); dist_pairs{3} = nan(m,m);
        for i = 1:m
            for j = 1:m
                for k = 1:m
                    dist_3(i,j,k) = compute_d_via_triangle_area([ subpoints( : , 1 , i ), subpoints( : , 2 , j ), subpoints( : , 3 , k ) ], eps);
                    dist_2(i,j,k) = compute_d_via_pairwise_sums([ subpoints( : , 1 , i ), subpoints( : , 2 , j ), subpoints( : , 3 , k ) ]);
                    dist_BaryC(i,j,k) = compute_d_via_barycenter([ subpoints( : , 1 , i ), subpoints( : , 2 , j ), subpoints( : , 3 , k ) ]);
                    dist_pairs{1}(i,j) = norm(subpoints( : , 1 , i ) -  subpoints( : , 2 , j ));
                    dist_pairs{2}(i,k) = norm(subpoints( : , 1 , i ) -  subpoints( : , 3 , k ));
                    dist_pairs{3}(j,k) = norm(subpoints( : , 2 , j ) -  subpoints( : , 3 , k ));
                end
            end
        end

        try
            switch (distance_func_type)
                case 1
                    probs = ones(3,m)/m;
                    [W, ~] = compute_W_n_dist_for_generic_dist_func(probs,dist_2);
                case 2
                    probs = ones(3,m)/m;
                    [W, ~] = compute_W_n_dist_for_generic_dist_func(probs,dist_3);
                case 3
                    probs = ones(3,m)/m;
                    [W, ~] = compute_W_n_dist_for_generic_dist_func(probs,dist_BaryC);
            end
        catch
            W = 10000;  % very rare, ignore hyperedge
        end

        if (isnan(W)) % very rare, ignore hyperedge
            W = 10000;
        end
        hyper_edges(t+sub_t-1,4) = W;
        
    end
        
    % now that we have the 4 hyperedges, we force them to violate
    % the triang ineq. % we perturbe as little as possible
    % we perturb the one hyperedge that is closest to violating the triang
    % ineq already
    
    a = hyper_edges(t,4); b = hyper_edges(t+1,4); c = hyper_edges(t+2,4); d = hyper_edges(t+3,4);
    [val , ix ] = min([(b+c+d) - a, (a+c+d) - b, (b+a+d) - c,(b+c+a) - d]);
    hyper_edges(t+ix-1,4) = hyper_edges(t+ix-1,4) + (1+strength_of_viol)*abs(val); % 30% violation %abs() not needed for n-metrics
    
end


fprintf("\n");

end



