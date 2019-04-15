options.R = 3;
options.d_ecc = 1;
options.h_ecc = 3;
options.theta_x = pi/6;
options.theta_y = pi/6;
options.theta_z = 0;

z = rot([0 1 0], options.theta_y)*rot([1 0 0], options.theta_x)*[0 0 1 1]';
zp = [z(1:2) ; 0];
phi = acos(zp' * [1 0 0]'/(norm(zp)+eps));
alpha = acos([0 0 1]*z(1:3)/(norm(z(1:3))+eps));

options.theta_x = 0;
options.theta_y = alpha;

figure
param = [0 0.2 5];
xx = param(1)*cos(phi) - param(2)*sin(phi);
yy = param(1)*sin(phi) + param(2)*cos(phi);
param(1:2) = [xx yy];

[flag, lowest_point, ~] = ARIE_findLowestPoint(param, options, 3);

X1 = X*cos(phi) - Y*sin(phi);
Y1 = X*sin(phi) + Y*cos(phi);