function mat_ix = convert_index_to_matrix(indices, num_classes)

    mat_ix = zeros(size(indices,1),num_classes);
    for i = 1:size(indices,1)
        mat_ix(i,indices(i)) = 1;
    end

end