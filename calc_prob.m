function [meshx, meshy, probs, prob] = calc_prob(points, input, avg, debug)
    a = 12;
    data_range = range([points, input], 2);
    data_mean = mean([points, input], 2);
    x_range = linspace(data_mean(1)-data_range(1), data_mean(1)+data_range(1), 200);
    y_range = linspace(data_mean(2)-data_range(2), data_mean(2)+data_range(2), 200);
    [meshx, meshy] = meshgrid(x_range, y_range);
    probs = zeros(size(meshx));
    prob = 0;
   
    for index = 1:length(points)
        x = points(1, index);
        y = points(2, index);
        if avg ~= 0
            rs = ((meshx - x).^2 + (meshy - y).^2).^0.5 / avg;
            C = log(a-1) / a;
            new_probs = (exp(-(rs + C - 1))./(exp(-a*(rs + C - 1)) + 1));
            r = ((input(1) - x).^2 + (input(2) - y).^2).^0.5 / avg;
            prob = prob + (exp(-(r + C - 1))./(exp(-a*(r + C - 1)) + 1));
        else
            rs = ((meshx - x).^2 + (meshy - y).^2).^0.5;
            new_probs = exp(-rs.^2/0.025);
            r = ((input(1) - x).^2 + (input(2) - y).^2).^0.5;
            prob = prob + exp(-r.^2/0.025);
        end
        probs = probs + new_probs;
    end
    
    if debug
        figure 
        hold on
        axis equal
        plot(points(1, :), points(2,:), 'bo')
        plot(input(1), input(2), 'rs')
        contour(meshx, meshy, probs);
        colorbar
        figure
        surf(meshx, meshy, probs);
        shading interp
        axis equal
    end
    
end
