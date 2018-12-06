% files = ["./Data/3-line-A.mat";...
%         "./Data/3-line-B.mat";...
%         "./Data/5-circle.mat";...
%         "./Data/5-circle-2.mat";...
%         "./Data/5-circle-3.mat";...
%         "./Data/5-circle-4.mat";...
%         "./Data/5-line-2.mat"];
    
% files = ["./Data/Course-1.mat"];
% files = ["./Data/3-line-A.mat";...
%         "./Data/3-line-B.mat";];
close all
xs = linspace(-5, 5, 100);
ys = linspace(-5, 5, 100);
rot = [1, 0; 0, 1];
trans = [0.05;0.05];

for index = 1:length(files)
    file = char(files(index, :));
    disp(file)
    new_scan = cluster_detection(file, 1);
    if index == 1
        cups = new_scan;
        avg_cup_dist = mean(norm(diff(cups, 1, 2)));
    else
        new_scan = generate_noise(new_scan, rot, trans);
        [cup_ind, corr_point_ind] = find_point_correspondence(cups, new_scan, avg_cup_dist);
        [R, t] = ICP(cups(:, cup_ind), new_scan(:, corr_point_ind));
        new_point_ind = setdiff(1:length(new_scan), corr_point_ind);
        new_points = new_scan(:, new_point_ind);
        cups = [cups, R * new_points + t];
    end
    z = calc_gradient(xs, ys, cups, 0.5);
    contour(xs, ys, z)
    
end

%%
function z = calc_gradient(xs, ys, points, order)
    [xmesh, ymesh] = meshgrid(xs, ys);
    z = zeros(size(xmesh));
    interpolated_points = 100;
    for i = 1:length(points) - 1
        x1 = points(1, i);
        y1 = points(2, i);
        x2 = points(1, i + 1);
        y2 = points(2, i + 1);
        plot([x1, x2], [y1, y2], 'k--')
        x_interp = linspace(x1, x2, interpolated_points);
        y_interp = linspace(y1, y2, interpolated_points);
        for j = 1:length(x_interp)
            z = z + (((xmesh - x_interp(j)).^2 + (ymesh - y_interp(j)).^2).^(-0.5 * order)) ./ interpolated_points;
        end
    end
end

function [R,t] = ICP(old_points, new_points)
    OP = old_points - mean(old_points, 2);
    figure
    hold on
    for interation = 1:3
        NP = new_points - mean(new_points, 2);
        covariance = NP * OP';
        [u, s, v] = svd(covariance);
        R = v*u';
        disp(R)
        t = mean(old_points, 2) - (R * mean(new_points, 2));
        disp(t)
        corrected_new_points = R*new_points + t
        all_points = [old_points, corrected_new_points];
        
        plot(old_points(1, :), old_points(2, :), 'bo')
        plot(new_points(1, :), new_points(2, :), 'rs')
        plot(corrected_new_points(1, :), corrected_new_points(2, :), '*')
        new_points = corrected_new_points;
    end
end

function points = generate_noise(input_points, rot, trans)
    noise = (rand(size(input_points))-0.5)*0.25;
    points = circshift(rot*(input_points+noise) + trans, 2);
end

function [corr_indices_global, corr_indices_local] = find_point_correspondence(global_points, encoder_corrected_points, avg_cup_dist)
    corr_indices_global = [];
    corr_indices_local = [];
    corr_dists = [];
    for index = 1:length(encoder_corrected_points)
        point = encoder_corrected_points(:, index);
        dists = vecnorm(global_points - point);
        [m, i] = min(dists);
        if m < avg_cup_dist * 0.4
            corr_indices_global = [corr_indices_global, i];
            corr_dists = [corr_dists, m];
            corr_indices_local = [corr_indices_local, index];
        end
    end
    [N, edges] = histcounts(corr_indices_global, length(global_points));
    dupes = find(N > 1);
    for global_dupe = dupes
        local_dupe_ind = find(corr_indices_global == global_dupe);
        [dist_to_keep, index_to_keep] = min(corr_dists(local_dupe_ind));
        local_ind_to_keep = corr_indices_local(index_to_keep);
        global_ind_to_keep = corr_indices_global(index_to_keep);
        corr_indices_local(local_dupe_ind) = [];
        corr_indices_local = [corr_indices_local, local_ind_to_keep];
        corr_indices_global(local_dupe_ind) = [];
        corr_indices_global = [corr_indices_global, global_ind_to_keep];
    end
end