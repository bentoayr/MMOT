function label_err = compute_labeling_error(clusters,num_classes,size_clust)

num_points = size(clusters,1);
label_err = inf;
for per = perms(1:num_classes)'
    M_per = eye(num_classes);
    M_per = M_per(1:num_classes,per);
    sol = kron(M_per,ones(size_clust,1));
    rr = sum(sum(abs(sol - clusters)));
    label_err = min(rr, label_err);
end
label_err = label_err/(2*num_points); 

end