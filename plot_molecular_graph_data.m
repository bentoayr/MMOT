%% make some plots of the data

tuncate_size = 20;


indices_bounds = [0, cumsum(num_per_data_set)];

hold on;
for i = 1:5
    
    Mat_form = zeros(length(indices_bounds(i)+1:indices_bounds(i+1)),tuncate_size);
        
    for j = indices_bounds(i)+1:indices_bounds(i+1)
        Mat_form(j - indices_bounds(i),:) = list_of_spects{j}(end-tuncate_size+1:end);
    end
    
    Mat_form = [real(Mat_form),imag(Mat_form)]; 
    
    [U , S, V] = svd(Mat_form);
    
    scatter3(U(:,1),U(:,2),U(:,3));
end
hold off;
