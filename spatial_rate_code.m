function [rate_matrix_total, x] = spatial_rate_code(matrix)

x = unique(matrix(~isnan(matrix(:,4)),4));

rate_matrix_total = [ ];

for i = 1 : length(x)
    [rate_matrix_smoothed_i, spike_count_i, spatial_occupancy_i] = rate_mtx5(matrix, x(i), 1, 30)
    rate_matrix_smoothed_i = rate_matrix_smoothed_i(:);
    rate_matrix_total(:,i) = [rate_matrix_smoothed_i];
   
end 

end 