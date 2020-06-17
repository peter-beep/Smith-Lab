function [rate_matrix, spike_count, spatial_occupancy] = rate_mtx5(dmtx, cluster, vis, bins)
% reduced data matrix for the particular visit
dmtx=dmtx(dmtx(:,6)==vis, :);

 
%all spike events in two vectors
xs = dmtx(dmtx(:,4)==cluster, 2);
ys = dmtx(dmtx(:,4)==cluster, 3);

%all time samples in two vectors 
 xt = dmtx(:,2);
 yt = dmtx(:,3);

%evenly spaced bins of x and y coordinate ranges
pos_edges_x = linspace(min(dmtx(:,2)), max(dmtx(:,2))+0.000000001, bins+1);
pos_edges_y = linspace(min(dmtx(:,3)), max(dmtx(:,3))+0.000000001, bins+1);


%2d histogram of event counts
spike_count = histcounts2(ys, xs, pos_edges_y, pos_edges_x);
spatial_occupancy = histcounts2(yt, xt, pos_edges_y, pos_edges_x)./100;

%flip
spike_count = flipud(spike_count);
spatial_occupancy = flipud(spatial_occupancy);

%divide spikes by time for rate
rate_matrix = spike_count./spatial_occupancy;