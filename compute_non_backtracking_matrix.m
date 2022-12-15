function [approx_non_back_track_matrix] = compute_non_backtracking_matrix(numV,adj)

    D = diag(sum(adj));
    approx_non_back_track_matrix = [adj, eye(numV) - D; eye(numV) , zeros(numV)];  
    
end