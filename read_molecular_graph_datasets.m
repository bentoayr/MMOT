function [list_of_spects,num_per_data_set ] =  read_molecular_graph_datasets(prune_flag,ignore_weights_flag,truncate_number)


dataset_names = ["dhfr_3d.sd" , "er_lit_3d.sd", "bzr_3d.sd", "cox2_3d.sd", "er_tox_3d.sd"];


count = 1;
list_of_spects = {};
num_per_data_set = [];
%graph_sizes = cell(length(dataset_names),1);


for dataset_name_id = 1:length(dataset_names)
    
    fid = fopen(strcat("./Datasets/",dataset_names(dataset_name_id)));
    
    count_per_data_set = 0;
    while(1)
        
        while(~feof(fid))
            res = fgetl(fid);
            if (res == "$$$$")
                break;
            end
        end
        if (feof(fid))
            break;
        end
        
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        %fgetl(fid);
        
        A = fscanf(fid, "%d ");
        numNodes = A(1);
        numEdges = A(2);
        for i = 1:numNodes+1
            fgetl(fid);
        end
        Adj_mat = zeros(numNodes);
        for i = 1:numEdges
            A = fscanf(fid, "%d",6);
            Adj_mat(A(1),A(2)) = A(3);
            Adj_mat(A(2),A(1)) = A(3);
        end
        
        
        if (ignore_weights_flag == 1)
            Adj_mat = 0.0 + (Adj_mat ~=0);
        end
        
        if (prune_flag == 1)
            node_deg = sum(Adj_mat~=0);
            non_leaf_nodes_ids = find(node_deg > 1);
            Adj_mat = Adj_mat(non_leaf_nodes_ids,non_leaf_nodes_ids);
            numNodes = size(Adj_mat,1);
            numEdges = sum(sum(Adj_mat~=0))/2;
        end
        
        if (prune_flag == 2)
            old_Adj_mat = 0;
            while( length(old_Adj_mat(:)) ~= length(Adj_mat(:)))
                old_Adj_mat = Adj_mat;
                node_deg = sum(Adj_mat~=0);
                non_leaf_nodes_ids = find(node_deg > 1);
                Adj_mat = Adj_mat(non_leaf_nodes_ids,non_leaf_nodes_ids);
                numNodes = size(Adj_mat,1);
                numEdges = sum(sum(Adj_mat~=0))/2;
            end
        end
        
        [approx_non_back_track_matrix] = compute_non_backtracking_matrix(numNodes,Adj_mat);
        spec_approx_non_back_track_matrix = sort(eig(approx_non_back_track_matrix)); %this spectrum differs from the previous one by +- 1 eigenvalues
        
        if (truncate_number >= 0)
            list_of_spects{count} = spec_approx_non_back_track_matrix(truncate_number+1:end);
            count = count + 1;
            count_per_data_set = count_per_data_set + 1;
        else
           
            if (numNodes == -truncate_number)
                list_of_spects{count} = spec_approx_non_back_track_matrix;
                count = count + 1;
                count_per_data_set = count_per_data_set + 1;
            end
            
        end
        
     
        %graph_sizes{dataset_name_id} = [graph_sizes{dataset_name_id}, numNodes];
            
    end

    
    num_per_data_set = [num_per_data_set, count_per_data_set];
end



