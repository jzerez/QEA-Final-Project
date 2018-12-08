%%Path plotting
clear all
close all

% test_cones = [5,-3;5,2;3,7;-2,7;-5,4;-5,-2;-3,-4;2,-4]';
% neato_origin = [4;-2;0];
test_cones = [1,1;2,1.5;3.75,3;4,5;4.5,7]';
neato_origin = [.5;1.25;0];
test_cones(3,:) = zeros(1,length(test_cones));
shuffled_cones = test_cones(:,randperm(length(test_cones)));
neato_orientation = [1;0;0];

figure
plot(test_cones(1,:),test_cones(2,:),'b-');
hold on
plot(shuffled_cones(1,:),shuffled_cones(2,:),'r-');
plot(neato_origin(1),neato_origin(2),'gs');
quiver(neato_origin(1),neato_origin(2),neato_orientation(1),neato_orientation(2),'AutoScale','off');
hold off
xlim([0,5])
axis equal
%% 

[~,closest_indx] = min(vecnorm(shuffled_cones-neato_origin))
ordered_cones(:,1) = shuffled_cones(:,closest_indx)
shuffled_cones = remove(ordered_cones(:,1),shuffled_cones);
ordered_cones(:,2) = next_cone(neato_origin,neato_orientation,shuffled_cones);
shuffled_cones = remove(ordered_cones(:,2),shuffled_cones);
ordered_cones(:,3) = next_cone(ordered_cones(:,2),ordered_cones(:,2)-ordered_cones(:,1),shuffled_cones);

next_pos = generate_offset(ordered_cones);
figure
plot(test_cones(1,:),test_cones(2,:),'b-');
hold on
plot(shuffled_cones(1,:),shuffled_cones(2,:),'r-');
plot(neato_origin(1),neato_origin(2),'gs');
quiver(neato_origin(1),neato_origin(2),neato_orientation(1),neato_orientation(2),'AutoScale','off');
plot(next_pos(1),next_pos(2),'ks');
quiver(neato_origin(1),neato_origin(2),next_pos(1)-neato_origin(1),next_pos(2)-neato_origin(2),'AutoScale','off');
hold off
xlim([0,5])
axis equal

k = cross(neato_orientation,next_pos-neato_origin);
angle3d = sign(k)*atan2d(norm(k),dot(neato_orientation,next_pos-neato_origin));
angle = angle3d(3);

forward = norm(next_pos-neato_origin);

% [func_points,func_paths] = generate_paths(test_cones);
% figure
% plot(test_cones(1,:),test_cones(2,:),'o');
% hold on
% plot(func_points(1,:),func_points(2,:),'s');
% quiver(func_points(1,1:length(func_paths)),func_points(2,1:length(func_paths)),func_paths(1,:),func_paths(2,:),'Color','r','LineStyle','--','AutoScale','off')
% hold off
% axis equal
%%
function [points,paths] = generate_paths(cones)
    vector = diff(cones,1,2);
    norm_vector = vector./vecnorm(vector);

    for i = 1:size(norm_vector,2)-1
        angle_between(i) = acosd(dot(norm_vector(:,i),norm_vector(:,i+1)));
        midangle(i) = 90+(angle_between(i)/2);
    end

    for j = 1:size(midangle,2)
        rot = [cosd(midangle(j)) -sind(midangle(j)) 0; sind(midangle(j)) cosd(midangle(j)) 0; 0 0 0];
        offset_dir(:,j) = rot*-0.25*norm_vector(:,j);
        points(1:2,j) = offset_dir(1:2,j)+cones(1:2,j+1);
    end
    paths = diff(points,1,2);
end

function point = generate_offset(cones)
    vector = diff(cones,1,2);
    norm_vector = vector./vecnorm(vector);

    for i = 1:size(norm_vector,2)-1
        angle_between(i) = acosd(dot(norm_vector(:,i),norm_vector(:,i+1)));
        midangle(i) = 90+(angle_between(i)/2);
    end

    for j = 1:size(midangle,2)
        rot = [cosd(midangle(j)) -sind(midangle(j)) 0; sind(midangle(j)) cosd(midangle(j)) 0; 0 0 0];
        offset_dir(:,j) = rot*-0.25*norm_vector(:,j);
        point(:,j) = offset_dir(:,j)+cones(:,j+1);
    end
end

function coneless = remove(cone,dataset)
    [~,indx] = max(sum(dataset==cone));
    dataset(:,indx) = []; 
    coneless = dataset;
end

function next = next_cone(origin,orientation,cones)
    vectors = cones-origin;
    distances = vecnorm(vectors);
    closest_indices = distances < 3;
    closest_cones = cones(:,closest_indices);
    closest_vectors = vectors(:,closest_indices);
    dot_length = 0;
    for i = 1:size(closest_vectors,2)
        if dot(closest_vectors(:,i),orientation) > dot_length
            dot_length = dot(closest_vectors(:,i),orientation);
            cone = closest_cones(:,i);
        end
    end
    next=cone;
end