% compute d via pairwise sums

function dist = compute_d_via_pairwise_sums(x)
    
    dist = norm(x(:,1) - x(:,2)) + norm(x(:,1) - x(:,3)) + norm(x(:,2) - x(:,3));
    
end