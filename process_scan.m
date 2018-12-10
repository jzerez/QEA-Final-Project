
% files = ["./Data/3-line-A.mat";...
%         "./Data/3-line-B.mat";...
%         "./Data/5-circle.mat";...
%         "./Data/5-circle-2.mat";...
%         "./Data/5-circle-3.mat";...
%         "./Data/5-circle-4.mat";...
%         "./Data/5-line-2.mat"];
%     
files = ["./Data/Course-1.mat";...
         "./Data/Course-2.mat";...
         "./Data/Course-3.mat"];
% files = ["./Data/3-line-A.mat";...
%         "./Data/3-line-B.mat";];
% close all
clf
xs = linspace(-5, 5, 100);
ys = linspace(-5, 5, 100);
lidar_offset = 0.1;         %m


sub = rossubscriber('/stable_scan');
lidar_pos = [0; 0];         %m, m
neato_pos = [lidar_offset; 0];
neato_ori = [1; 0];         %m, m

hold on
new_scan = cluster_detection(take_scan(sub), 1);
pause
cups = new_scan;
avg_cup_dist = mean(norm(diff(cups, 1, 2)));
plot_elems = [];
 for index = 1:55
     % Store last position info
    last_lidar_pos = lidar_pos;
    last_neato_pos = neato_pos;
    
    % Move
    [angle, translation] = calcmove(neato_pos, neato_ori, cups); 
     pause(0.5)
    
    % Update position info
    neato_pos = last_neato_pos + translation;
    neato_ori = (neato_pos - last_neato_pos) / norm(neato_pos - last_neato_pos);
    lidar_pos = neato_pos - (lidar_offset * neato_ori);
    
    % Scan
    new_scan = cluster_detection(take_scan(sub), 0);
    
    % Plot
    plot_elems(1) = plot(lidar_pos(1), lidar_pos(2), 'go');
    quiver(lidar_pos(1), lidar_pos(2), neato_ori(1)*0.2, neato_ori(2)*0.2, 'g')
    plot_elems(2) = plot(cups(1, :), cups(2, :), 'r*');
    
    % Translate based on encoders
    new_scan = rough_transform(neato_ori, lidar_pos, new_scan);
    % Bayesian correspondence
    [cup_ind, corr_point_ind] = calc_bayes_corr(cups, new_scan, avg_cup_dist);
    plot_elems(3) = plot(new_scan(1,:), new_scan(2, :), 'g*');
    plot_elems(4) = plot(new_scan(1, corr_point_ind), new_scan(2, corr_point_ind), 'ks');
    % ICP
    [R, t] = ICP(cups(:, cup_ind), new_scan(:, corr_point_ind));
    new_scan = R * new_scan + t;
    lidar_pos = R * lidar_pos + t;
    
    [cup_ind, corr_point_ind] = calc_bayes_corr(cups, new_scan, avg_cup_dist);
    [R, t] = ICP(cups(:, cup_ind), new_scan(:, corr_point_ind));
    new_scan = R * new_scan + t;
    lidar_pos = R*lidar_pos + t;
    neato_ori = (lidar_pos - last_lidar_pos) / norm(lidar_pos - last_lidar_pos);
    
    [cup_ind, corr_point_ind] = calc_bayes_corr(cups, new_scan, avg_cup_dist);
    [R, t] = ICP(cups(:, cup_ind), new_scan(:, corr_point_ind));
  
    cups(:, cup_ind) = (cups(:, cup_ind) + new_scan(:, corr_point_ind)) / 2;
    
    new_point_ind = setdiff(1:length(new_scan), corr_point_ind);
    new_points = new_scan(:, new_point_ind);
    cups = [cups, R * new_points + t];
    avg_cup_dist = mean(norm(diff(cups, 1, 2)));
    plot_elems(5) = plot(cups(1,:), cups(2, :), 'k*');

    lidar_pos = R*lidar_pos + t;
    neato_pos = lidar_pos + (lidar_offset * neato_ori);
    neato_ori = (lidar_pos - last_lidar_pos) / norm(lidar_pos - last_lidar_pos);
    plot_elems(6) = plot(lidar_pos(1), lidar_pos(2), 'bo');
    quiver(lidar_pos(1), lidar_pos(2), neato_ori(1)*0.2, neato_ori(2)*0.2, 'b')
    legend(plot_elems, {'encoder pos', 'old cones' 'new cones', 'corr cones', 'updated cones', 'final lidar'})
    pause
    clf
    hold on
    axis equal
    plot(cups(1,:), cups(2, :), 'k*');
    plot(lidar_pos(1), lidar_pos(2), 'ro');
    plot(neato_pos(1), neato_pos(2), 'rs');
end


%%
function polar_points = take_scan(sub)
    scan_message = receive(sub);
    r = scan_message.Ranges(1:end-1);
    theta = [0:359]';
    polar_points = [r, theta];
end
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

function [R_final,T_final] = ICP(old_points, new_points)
    OP = old_points - mean(old_points, 2);
    R_final = eye(2);
    T_final = [0;0];
    for interation = 1:3
        NP = new_points - mean(new_points, 2);
        covariance = NP * OP';
        [u, ~, v] = svd(covariance);
        reflect = eye(length(v));
        reflect(end) = det(v*u');
        R = v*reflect*u';
        disp(R)
        t = mean(old_points, 2) - (R * mean(new_points, 2));
        disp(t)
        corrected_new_points = R*new_points + t;
        all_points = [old_points, corrected_new_points];
        
        new_points = corrected_new_points;
        R_final = R_final * R;
        T_final = T_final + t;
    end
end

function points = generate_noise(input_points, rot, trans)
    noise = (rand(size(input_points))-0.5)*0.25;
    points = circshift(rot*(input_points+noise) + trans, 2);
end

function points = rough_transform(orientation, trans, input_points)
    angle = atan2d(orientation(2), orientation(1));
    rot = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
    points = rot * input_points + trans;
end

function [pos, ori, angle] = update_pos_and_ori(angle, translation, old_pos, old_ori, old_angle)
    rot = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
    pos = old_pos + translation;
    ori = rot * old_ori;
    angle = angle + old_angle;
    
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

function [corr_indices_global, corr_indices_local] = calc_bayes_corr(global_points, encoder_corrected_points, avg_cup_dist)
    pr_new = 0.25;
    pr_old = 1-pr_new;
    corr_indices_local = [];
    corr_indices_global = [];
    for index = 1:length(encoder_corrected_points)
        point = encoder_corrected_points(:, index);
        [~, ~, new_probs, new_prob] = calc_prob(global_points, point, avg_cup_dist, 0);
        [X, Y, old_probs, old_prob] = calc_prob(global_points, point, 0, 0);
        unnormalized_new_probability = pr_new * new_prob;
        unnormalized_old_probability = pr_old * old_prob;
        posterior_new = unnormalized_new_probability / (unnormalized_new_probability + unnormalized_old_probability);
        if posterior_new < 0.5
            corr_indices_local = [corr_indices_local; index];
            [M, I] = min(vecnorm(global_points - point));
            corr_indices_global = [corr_indices_global; I];
        end
%         contour(X, Y, old_probs)
    end
end

