function [clusters_nhcut,clusters_ttm] = Hypergraph_clustering_th(num_points,hyper_edges,threshold,max_num_clust)

clusters_nhcut = randi(max_num_clust,num_points,1); clusters_nhcut = convert_index_to_matrix(clusters_nhcut, max_num_clust);
clusters_ttm = randi(max_num_clust,num_points,1); clusters_ttm = convert_index_to_matrix(clusters_ttm, max_num_clust);

num_hyp_edges = size(hyper_edges,1);
H = inf(num_points,num_hyp_edges);

for i = 1:num_hyp_edges
    H(hyper_edges(i,1:3),i) = hyper_edges(i,4);
end

N = num_points;
edgwt = hyper_edges(:,4);
H = exp(-H) > threshold;
K = max_num_clust;

try
    W = (H.*repmat(edgwt',N,1))*(H'./3);
    W(isnan(W)) = 0;
    D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
    W = eye(N) - W.*dd;
    [vec,val] = eig(W,'nobalance');
    temp = sortrows([diag(val) vec'],1);
    evecs = temp(1:K,2:end)';
    for i = 1:N
        if (norm(evecs(i,:))>0)
            evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
        end
    end
    clusters_nhcut = kmeans(evecs,K,'emptyaction','singleton','replicates',30);
    clusters_nhcut = convert_index_to_matrix(clusters_nhcut, K);
catch
end

A = zeros(N,N,N);
for i = 1:num_hyp_edges
    if (  exp(-hyper_edges(i,4 ) ) > threshold)
        i1 = hyper_edges(i,1);
        i2 = hyper_edges(i,2);
        i3 = hyper_edges(i,3);
        
        A(i1,i2,i3) = 1; A(i1,i3,i2) = 1;
        A(i2,i1,i3) = 1; A(i2,i3,i1) = 1;
        A(i3,i2,i1) = 1; A(i3,i1,i2) = 1;
    end
end

try
    W = sum(A,3);
    W(isnan(W)) = 0;
    D = sqrt(sum(W,2)); D(D==0) = 1; dd = 1./(D*D');
    W = W.*dd;
    [vec,val] = eig(W,'nobalance');
    temp = sortrows([diag(val) vec'],-1);
    evecs = temp(1:K,2:end)';
    for i = 1:N
        if (norm(evecs(i,:))>0)
            evecs(i,:) = evecs(i,:)./norm(evecs(i,:));
        end
    end
    clusters_ttm = kmeans(evecs,K,'emptyaction','singleton','replicates',30);
    clusters_ttm = convert_index_to_matrix(clusters_ttm, K);
catch
end

end





