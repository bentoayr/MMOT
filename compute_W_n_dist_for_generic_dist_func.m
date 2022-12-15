function [W,p] = compute_W_n_dist_for_generic_dist_func(probs,dist)

    m = size(probs,2); %number of points in each distribution
 
    cvx_begin quiet
        variable p(m,m,m)
        
        minimize    (  sum(sum(sum(p.*(dist)))  ) )
        subject to
            p >= 0;
            sum(sum(sum(p))) == 1;
            sum(sum(permute(p,[1,2,3]),3),2) == probs(1,:)';
            sum(sum(permute(p,[2,1,3]),3),2) == probs(2,:)';
            sum(sum(permute(p,[3,1,2]),3),2) == probs(3,:)';
    
    cvx_end
    
    if (isnan(cvx_optval) || isinf(abs(cvx_optval)))  % if CVX did not work, try a small perturbation.
        cvx_begin quiet
            variable p(m,m,m)

            minimize    (  sum(sum(sum(p.*(dist+ 0.05*max(dist(:))*rand(m,m,m))))  ) )
            subject to
            p >= 0;
            sum(sum(sum(p))) == 1;
            sum(sum(permute(p,[1,2,3]),3),2) == probs(1,:)';
            sum(sum(permute(p,[2,1,3]),3),2) == probs(2,:)';
            sum(sum(permute(p,[3,1,2]),3),2) == probs(3,:)';
        
        cvx_end
        
    end
    
    W = cvx_optval;

end