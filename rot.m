function T = rot(f, theta)
% Perform general rotation transformation
%
% Parameters:
%  - f: rotation vector
%  - theta: rotation angle (radian)
%
% Author:
%  Rui Li (raysworld@outlook.com)
%
% Date:
%  2014-06-14 File created
%  2014-07-17 Add comments
%

% set coefficients
kx = f(1);
ky = f(2);
kz = f(3);
V = 1 - cos(theta);

% set the g.r.t matrix
T = zeros(4, 4);
T(1,1) = kx*kx*V + cos(theta);
T(1,2) = ky*kx*V - kz*sin(theta);
T(1,3) = kz*kx*V + ky*sin(theta);
T(2,1) = kx*ky*V + kz*sin(theta);
T(2,2) = ky*ky*V + cos(theta);
T(2,3) = kz*ky*V - kx*sin(theta);
T(3,1) = kx*kz*V - ky*sin(theta);
T(3,2) = ky*kz*V + kx*sin(theta);
T(3,3) = kz*kz*V + cos(theta);
T(4,4) = 1;
