files = ["./Data/5-circle.mat";...
        "./Data/5-circle-2.mat";...
        "./Data/5-circle-3.mat";...
        "./Data/5-circle-4.mat";...
        "./Data/5-line-2.mat"];
% files = ["./Data/5-line-2.mat"];
close all

for index = 1:length(files)
    file = char(files(index, :));
    disp(file)
    cups = cluster_detection(file, 1);
    xs = linspace(-5, 5, 100);
    ys = linspace(-5, 5, 100);
    z = calc_gradient(xs, ys, cups, 0.5);
    contour(xs, ys, z)
    
    rot = [0, 1; -1, 0];
    trans = [2;3];
    noise = (rand(size(cups))-0.5)*0.25;
    new_points = circshift(rot*(cups+noise) + trans, 2);
    ICP(cups, new_points)
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
        x_interp = linspace(x1, x2, interpolated_points);
        y_interp = linspace(y1, y2, interpolated_points);
        for j = 1:length(x_interp)
            z = z + (((xmesh - x_interp(j)).^2 + (ymesh - y_interp(j)).^2).^(-0.5 * order)) ./ interpolated_points;
        end
    end
end

function all_points = ICP(old_points, new_points)
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