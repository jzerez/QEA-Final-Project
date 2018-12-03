%Path plotting
%assuming counter clockwise motion
clear all
close all

test_cones = [5,-3;5,2;3,7;-2,7;-5,4;-5,-2;-3,-4;2,-4;5,-3];
z = zeros(length(test_cones),1);
test_cones_3d = cat(2,test_cones,z);
figure
plot(test_cones(:,1),test_cones(:,2),'o');

for i = 1:length(test_cones_3d)-1
    vector(i,:) = test_cones_3d(i+1,:)-test_cones_3d(i,:);
end

for j = 1:length(vector)-1
    angle_between(j) = atan2d(norm(cross(vector(j,:),vector(j+1,:))),dot(vector(j,:),vector(j+1,:)));
    midangle(j) = (180-angle(j))/2;
    midangle(j) = (atan2d(norm(cross(vector(j,:),vector(j+1,:))),dot(vector(j,:),vector(j+1,:)))+180)/2;
end

