files = ["./Data/Course-1.mat";...
         "./Data/Course-2.mat";...
         "./Data/Course-3.mat"];

clf
xs = linspace(-5, 5, 100);
ys = linspace(-5, 5, 100);
lidar_offset = 0.1;         %m

sub = rossubscriber('/stable_scan');
neato_pos = [0; 0];         %m, m
neato_ori = [1; 0];         %m, m

hold on
new_scan = cluster_detection(take_scan(sub), 1) - (lidar_offset * neato_ori);
pause
cups = new_scan;
avg_cup_dist = calc_avg_dist(cups);
plot_elems = [];
 for index = 1:55
     % Store last position info
    last_neato_pos = neato_pos;
    
    % Move
    [angle, translation] = calcmove(neato_pos, neato_ori, cups); 
     pause(1)
    
    % Update position info
    neato_pos = last_neato_pos + translation;
    neato_ori = (neato_pos - last_neato_pos) / norm(neato_pos - last_neato_pos);
    
    % Scan
    new_scan = cluster_detection(take_scan(sub), 0) - (lidar_offset * neato_ori);
    
    % Plot
    plot_elems(1) = plot(neato_pos(1), neato_pos(2), 'go');
    quiver(neato_pos(1), neato_pos(2), neato_ori(1)*0.2, neato_ori(2)*0.2, 'g')
    plot_elems(2) = plot(cups(1, :), cups(2, :), 'r*');
    
    % Translate based on encoders
    new_scan = rough_transform(neato_ori, neato_pos, new_scan);
    % Bayesian correspondence
    [cup_ind, corr_point_ind] = calc_bayes_corr(cups, new_scan, avg_cup_dist, neato_pos);
    plot_elems(3) = plot(new_scan(1,:), new_scan(2, :), 'g*');
    plot_elems(4) = plot(new_scan(1, corr_point_ind), new_scan(2, corr_point_ind), 'ks');
    % ICP
    [R, t] = ICP(cups(:, cup_ind), new_scan(:, corr_point_ind));
    new_scan = R * new_scan + t;
    neato_pos = R * neato_pos + t;
    plot(neato_pos(1), neato_pos(2), 'b+')
    % Do correspondence and ICP two more times
    for jindex = 1:2
        [cup_ind, corr_point_ind] = calc_bayes_corr(cups, new_scan, avg_cup_dist, neato_pos);
        [R, t] = ICP(cups(:, cup_ind), new_scan(:, corr_point_ind));
        new_scan = R * new_scan + t;
         neato_pos = R*neato_pos + t;
        plot(neato_pos(1), neato_pos(2), 'b+')
    end
    neato_ori = (neato_pos - last_neato_pos) / norm(neato_pos - last_neato_pos);  
    % Average corresponding cups
    cups(:, cup_ind) = (cups(:, cup_ind) + new_scan(:, corr_point_ind)) / 2;
    
    % Add new cups
    new_point_ind = setdiff(1:length(new_scan), corr_point_ind);
    new_points = new_scan(:, new_point_ind);
    cups = [cups, new_points];
    calc_avg_dist(cups);
    
    % Plot
    plot_elems(5) = plot(cups(1,:), cups(2, :), 'k*');
    plot_elems(6) = plot(neato_pos(1), neato_pos(2), 'bo');
    quiver(neato_pos(1), neato_pos(2), neato_ori(1)*0.2, neato_ori(2)*0.2, 'b')
    legend(plot_elems, {'enco    der pos', 'old cones' 'new cones', 'corr cones', 'updated cones', 'final lidar'})
    pause
    clf
    hold on
    axis equal
    plot(cups(1,:), cups(2, :), 'k*');
    plot(neato_pos(1), neato_pos(2), 'ro');
    plot(neato_pos(1), neato_pos(2), 'rs');
end


%%  
function polar_points = take_scan(sub)
    scan_message = receive(sub);
    r = scan_message.Ranges(1:end-1);
    theta = [0:359]';
    polar_points = [r, theta];
end

function avg = calc_avg_dist(points)
    min_dists = zeros([1, size(points, 2)]);
    for index = 1:size(points, 2)
        point = points(:, index);
        unique = points;
        unique(:, index) = [];
        dists = vecnorm(unique - point);
        min_dists(index) = min(dists);
    end
    avg = mean(min_dists);
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

function [corr_indices_global, corr_indices_local] = calc_bayes_corr(global_points, encoder_corrected_points, avg_cup_dist, neato_pos)
    pr_new = 0.2;
    pr_old = 1-pr_new;
    corr_indices_local = [];
    corr_indices_global = [];
    close_indices = find(vecnorm(global_points - neato_pos) < 1.5);
    close_global_points = global_points(:,  close_indices); 
    viscircles(close_global_points', 0.03*ones(1, size(close_global_points, 2)))
    for index = 1:length(encoder_corrected_points)
        point = encoder_corrected_points(:, index);
        [~, ~, new_probs, new_prob] = calc_prob(close_global_points, point, avg_cup_dist, 0);
        [X, Y, old_probs, old_prob] = calc_prob(close_global_points, point, 0, 0);
        unnormalized_new_probability = pr_new * new_prob;
        unnormalized_old_probability = pr_old * old_prob;
        posterior_new = unnormalized_new_probability / (unnormalized_new_probability + unnormalized_old_probability);
        if posterior_new < 0.6 
            corr_indices_local = [corr_indices_local; index];
            [~, I] = min(vecnorm(close_global_points - point));
            corr_indices_global = [corr_indices_global; close_indices(I)];
        else
            disp("NEW POINT: ")
            disp(index)
        end
        contour(X, Y, new_probs, '--' )
        contour(X, Y, old_probs)
    end
end