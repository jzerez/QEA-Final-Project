%rosinit('10.0.75.2',11311, 'NodeHost','10.0.75.1')
clear all
sub = rossubscriber('/stable_scan');

% Collect data at the room origin
scan_message = receive(sub);
r_1 = scan_message.Ranges(1:end-1);
theta_1 = [0:359]';
scan1 = [r_1, theta_1];
figure
polarplot(deg2rad(theta_1), r_1, 'bo')
