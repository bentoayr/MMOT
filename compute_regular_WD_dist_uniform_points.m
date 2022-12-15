% this uses the norm between points as distances
function [W,p] = compute_regular_WD_dist_uniform_points(distrib1,distrib2)

    n1 = size(distrib1,2);
    n2 = size(distrib2,2);
    
    dist = nan(n1,n2);
    for i  = 1:n1
        for j = 1:n2
            dist(i,j) = norm( distrib1(:,i) - distrib2(:,j));
        end
    end

    probs{1} = ones(n1,1)/n1;
    probs{2} = ones(n2,1)/n2;

    cvx_begin quiet

    variable ppair(n1,n2)

    minimize    (  sum(sum(ppair.*dist))     )  
    subject to
        ppair >= 0;
        sum(sum(ppair)) == 1;
        sum(ppair,1) == probs{2}';
        sum(ppair,2) == probs{1};

    cvx_end

    if (isnan(cvx_optval) || isinf(abs(cvx_optval)) ) % if CVX did not work, try a small perturbation.
        cvx_begin quiet
        variable ppair(n1,n2)
        minimize    (  sum(sum(ppair.*(dist + 0.05*max(dist(:))*rand(n1,n2)   )))     )  
        subject to
            ppair >= 0;
            sum(sum(ppair)) == 1;
            sum(ppair,1) == probs{2}';
            sum(ppair,2) == probs{1};

        cvx_end
   
    end
    
    
    if (isnan(cvx_optval) || isinf(abs(cvx_optval)) ) % if for some reason cvx fails, output a large distance
        W = 1000000;
    else
        W = cvx_optval;
    end
    
    p = ppair;
    

end