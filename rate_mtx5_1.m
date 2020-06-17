function [rate_matrix_smoothed, spike_count, spatial_occupancy] = rate_mtx5_1(dmtx, cluster, vis, bins)
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

% To get rid of the interneurons with an average firing rate higher than 5Hz
TimeLength = max(dmtx(:,1))-min(dmtx(:,1));
MeanFreq = sum(spike_count(:))/TimeLength;

%divide spikes by time for rate
if ~ismember(cluster, unique(dmtx(:,4)))
rate_matrix_smoothed = spike_count./spatial_occupancy;

end



if ismember(cluster, unique(dmtx(:,4)))

% if MeanFreq < 5
rate_matrix_smoothed = skagg_smooth(spike_count, spatial_occupancy);
%rate_matrix_smoothed = smooth2a(rate_matrix, size(rate_matrix,1)/10, size(rate_matrix,2)/10);

rate_matrix_smoothed(isinf(rate_matrix_smoothed)) = 0 ;
 else
   rate_matrix_smoothed = nan(20,20);
end
%INTERNAL FUNCTION
%
%adaptive smoothing from Skaggs, McNaughton, Wilson, & Barnes, 1996
    function mtx_smth = skagg_smooth(spike_counts, occupancy_counts)
    % An adaptive smoothing method to optimize the trade-off between blurring 
    % error and sampling error. The firing rate at each bin was estimated
    % by expanding a circle around the point until the radius of the circle 
    % (in bins) is greater than a constant (try 10000) divided by n*sqrt(s),
    % where n is the number of occupancy samples within the circle and s is 
    % the total number of spikes within the circle. With a position sampling 
    % rate of 100 Hz, the firing rate at that point was then set to 100*n*s.
    %
    %Skaggs, McNaughton, Wilson, Barnes 1996, see also Ito 2015. This method is
    %favored by the Mosers.

    %skagg constant
    skg_c = 10000;

    occupancy_counts = occupancy_counts.*100;
    occupancy_counts(occupancy_counts==0) = nan;

    %preallocate rate matrix
    mtx_smth = nan(size(spike_counts));

        %iterate through rows
        for ir = 1:size(spike_counts,1)

            %iterate through columns
            for ic = 1:size(spike_counts,2)

                %skip (leave nan) if no occupancy
                if isnan(occupancy_counts(ir, ic))
                    continue
                end

                %preset test values counter
                radius = 0;
                skagg_val = 1; %arbitrarily higher than radius

                %keep trying until pass skagg test
                while radius < skagg_val

                    %expand smooth mask size
                    radius = radius+1;

                    %Circle with radius centered at ir,ic
                    [mgx, mgy] = meshgrid(1:size(spike_counts,1), 1:size(spike_counts,2));
                    circle_idx = sqrt((mgx-ic).^2+(mgy-ir).^2)<=radius;

                    %sum within circle area
                    circ_spikes = nansum(spike_counts(circle_idx));
                    circ_occupancy = nansum(occupancy_counts(circle_idx));

                    %calculate skagg test
                    skagg_val = skg_c / (circ_occupancy * sqrt(circ_spikes));

                end

                %set smoothed rate
                mtx_smth(ir, ic) = 100 * circ_spikes / circ_occupancy;
                
%                ir
%                 ic
%                 circ_spikes
%                 circ_occupancy
%                 100 * circ_spikes / circ_occupancy
                
            end

        end
    end
end