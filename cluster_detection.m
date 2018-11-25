close all
clear all
file = 'Data/5-circle.mat';
load(file)

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

polarplot(thetas, rs, 'bo')
figure
plot(delta_rs)

clusters = [];
first_index = 1;
for index = 1:length(delta_rs)
    if abs(delta_rs(index)) > radius_std
        clusters = [clusters, [first_index; index]];
        first_index = index + 1; 
    end
end

[xs, ys] = pol2cart(thetas, rs);
cluster_centers = zeros(size(clusters));
for index = 1:length(clusters)
   first = clusters(1, index);
   last = clusters(2, index);
   xavg = mean(xs(first:last));
   yavg = mean(ys(first:last));
   cluster_centers(:, index) = [xavg; yavg];
end
figure
hold on
plot(xs, ys, 'bo')
plot(0, 0, 'gs')
axis equal
for center = 1:length(clusters)
    diff(clusters(:, center))
    plot(cluster_centers(1, center), cluster_centers(2, center), 'k*')
end

cup_indicies = [];
for index = 1:length(clusters)
    r = norm(cluster_centers(:, index));
    theta = (abs(diff(clusters(:, index))) + 1) / 360 * 2 * pi;
    clusters(:, index)
    s = r * theta;
    if s > DESIRED_RADIUS * 0.5 && s < DESIRED_RADIUS * 2
        plot(cluster_centers(1, index), cluster_centers(2, index), 'rs')
        cup_indicies = [cup_indicies, index]; 
    end
end
