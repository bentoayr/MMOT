% compute d via barycenter 

function dist = compute_d_via_barycenter(x)
    
    p.Data = x';
    med = Weiszfeld(p);
    
    dist = norm(x(:,1) - med.xMedian) + norm(x(:,2) - med.xMedian) + norm(x(:,3) - med.xMedian);
    
    if (isnan(dist)) % if rare problem with Weiszfeld, return approximate heuristic
        med = mean(x,2); 
        dist = norm(x(:,1) - med) + norm(x(:,2) - med) + norm(x(:,3) - med);
    end
    
end