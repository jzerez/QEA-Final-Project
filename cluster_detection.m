function cones = cluster_detection(file, debug)
    disp(file)
    load(char(file))

    DESIRED_RADIUS = 0.0381;   %m

    rs = scan1(:, 1);
    thetas = deg2rad(scan1(:, 2));
    non_zero = find(rs ~= 0);
    rs = rs(non_zero);
    thetas = thetas(non_zero);
    [thetas, sorted_inds] = sort(thetas);
    rs = rs(sorted_inds);

    delta_rs = diff(rs);
    radius_std = std(delta_rs) * 0.2;
    clusters = [];
    first_index = 0;
    for index = 1:length(delta_rs)
        if abs(delta_rs(index)) > radius_std
            rs = circshift(rs, index);
            thetas = circshift(thetas, index);
            delta_rs = diff(rs);
            first_index = 1;
            break
        end
    end
    for index = 1:length(delta_rs)
        if abs(delta_rs(index)) > radius_std
            clusters = [clusters, [first_index; index]];
            first_index = index + 1;
        end
    end
    % for index = 1:length(delta_rs) - 1
    %     secondary_index = mod(first_index - index + length(delta_rs), length(delta_rs)); ;
    %     if abs(delta_rs(secondary_index)) > radius_std
    %         clusters = [clusters, [first_index; secondary_index]];
    %         break
    %     end
    % end

    [xs, ys] = pol2cart(thetas, rs);
    cluster_centers = zeros(size(clusters));
    for index = 1:length(clusters)
       first = clusters(1, index);
       last = clusters(2, index);
       xavg = mean(xs(first:last));
       yavg = mean(ys(first:last));
       cluster_centers(:, index) = [xavg; yavg];
    end

    
    cup_indicies = [];
    for index = 1:length(clusters)
        r = norm(cluster_centers(:, index));
        theta = (abs(diff(clusters(:, index))) + 1) / 360 * 2 * pi;
        clusters(:, index)
        s = r * theta;
        if s > DESIRED_RADIUS * 0.75 && s < DESIRED_RADIUS * 2
            cup_indicies = [cup_indicies, index]; 
        end
    end
    cones = [cluster_centers(1, cup_indicies); cluster_centers(2, cup_indicies)];


     if debug
        figure
        hold on
        plot(xs, ys, 'bo')
        plot(0, 0, 'gs')
        axis equal
        for center = 1:length(clusters)
            diff(clusters(:, center))
            plot(cluster_centers(1, center), cluster_centers(2, center), 'k*')
        end
        
        plot(cones(1, :), cones(2, :), 'r-') 
     end
end