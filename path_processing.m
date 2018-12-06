%Path plotting
%assuming counter clockwise motion
clear all
close all

test_cones = [5,-3;5,2;3,7;-2,7;-5,4;-5,-2;-3,-4;2,-4;5,-3];

for i = 1:length(test_cones)-1
    vector = diff(test_cones)';
    norm_vector = vector./vecnorm(vector);
end

for j = 1:length(norm_vector)-1
    angle_between(j) = acosd(dot(norm_vector(:,j),norm_vector(:,j+1)));
    midangle(j) = 90+(angle_between(j)/2);
end

for k = 1:length(midangle)
    rot = [cosd(midangle(k)) -sind(midangle(k)); sind(midangle(k)) cosd(midangle(k))];
    offset_dir(:,k) = rot*norm_vector(:,k);
    offset_point(:,k) = offset_dir(:,k)+test_cones(k+1,:)';
end

for m = 1:length(offset_point)-1
    path_vector(:,m) = offset_point(:,m+1)-offset_point(:,m);
end

figure
plot(test_cones(:,1),test_cones(:,2),'o');
hold on
plot(offset_point(1,:),offset_point(2,:),'s');
quiver(test_cones(1:8,1)',test_cones(1:8,2)',vector(1,:),vector(2,:),'Color','b','LineStyle','--','AutoScale','off')
quiver(offset_point(1,1:6),offset_point(2,1:6),path_vector(1,:),path_vector(2,:),'Color','r','LineStyle','--','AutoScale','off')
hold off
axis equal
