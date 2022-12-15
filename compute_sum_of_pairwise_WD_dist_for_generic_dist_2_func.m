function [W,p] = compute_sum_of_pairwise_WD_dist_for_generic_dist_2_func(probs,dist)

 
    W = 0;
    for i = 1:length(dist) %sum over all pairs
       
        dist_pair = dist{i};
        probs_pair = probs{i};
        
        m = size(probs_pair,2); %number of points in each distribution
        
        cvx_begin quiet
    
        variable ppair(m,m)
        
        minimize    (  sum(sum(ppair.*dist_pair))     )  
        subject to
            ppair >= 0;
            sum(sum(ppair)) == 1;
            sum(ppair,1) == probs_pair(2,:);
            sum(ppair,2) == probs_pair(1,:)';
            
        cvx_end
    
        W = W + cvx_optval;
        p{i} = ppair;
    end
    
    

end