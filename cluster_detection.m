function cone_centers = cluster_detection(scan1, debug)
    DESIRED_RADIUS = 0.0351;   %m
    rs = scan1(:, 1);
    thetas = deg2rad(scan1(:, 2));
    [thetas, rs] = cleanData(thetas, rs);
   
    delta_rs = diff(rs);
    radius_std = std(delta_rs) * 0.2;
    clusters = [];
    first_index = 0;
    for index = 1:length(delta_rs)
        if abs(delta_rs(index)) > radius_std
            if first_index == 0
                first_index = index + 1;
            else
                clusters = [clusters, [first_index; index]];
                first_index = index + 1;
            end
        end
    end
    rs = circshift(rs, -clusters(1,1) + 1);
    thetas = circshift(thetas, -clusters(1,1)+1);
    delta_rs = diff(rs);
    clusters = clusters - [clusters(1,1)-1; clusters(1,1)-1];
    
    for index = clusters(end)+1:length(delta_rs)
        if abs(delta_rs(index)) > radius_std
            clusters = [clusters, [clusters(end)+1; index]];
        end
    end
    clusters = [clusters, [clusters(end)+1; length(rs)]];
    
    large_clusters = diff(clusters) > 2;
    clusters = clusters(:, large_clusters);
    cones_indices = [];
    cone_centers = [];
    [xs, ys] = pol2cart(thetas, rs);
    
    for index = 1:size(clusters, 2)
        start_index = (clusters(1, index));
        end_index = (clusters(2, index));
        sample_points = round(linspace(start_index, end_index, 3));
        cluster_points = [xs(sample_points)'; ys(sample_points)'];
        [center, radius] = calc_circle(cluster_points);
        if debug
            plot(xs(start_index:end_index), ys(start_index:end_index), '*')
        end
        avg_radius = mean(vecnorm(cluster_points - center));
        if avg_radius > DESIRED_RADIUS * 0.4 && avg_radius < DESIRED_RADIUS * 2
            cones_indices = [cones_indices, index]; 
            cone_centers = [cone_centers, center];
        end
    end
    
    if debug
        plot(0, 0, 'gs')
        if numel(cone_centers) > 0
            plot(cone_centers(1, :), cone_centers(2, :), 'r*') 
        else
            disp("WARNING: NO CONES")
        end
        axis equal
    end


    
     
    function [center,r] = calc_circle(points)
        p1 = points(:, 1);
        p2 = points(:, 2);
        p3 = points(:, 3);
        mid1 = (p2+p1) / 2;
        mid2 = (p3+p2) / 2;
        slope1 = -(p2(1) - p1(1)) / (p2(2)-p1(2));
        slope2 = -(p3(1) - p2(1)) / (p3(2)-p2(2));
        system = rref([1, -slope1, mid1(2) - (slope1 * mid1(1));
                  1, -slope2, mid2(2) - (slope2 * mid2(1))]);
        center = [system(2, 3);system(1,3)];
        r = norm(center-p1);
    end

    function [ctheta,cr] = cleanData(theta, r)
        nonzero_r = r ~= 0;
        close_r = r < 1.5;
        i_clean = nonzero_r & close_r;  % indices of clean data
        ctheta = theta(i_clean);
        cr = r(i_clean);
    end
     
end