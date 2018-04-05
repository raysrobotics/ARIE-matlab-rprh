% This .m script is used for visualizing the ARIE 
% in the assemby task of round-peg and round-hole
%
% This version (script_main_rot) fixed the problem (in some way)  
%   - when theta_x ~= 0 and theta_y ~= 0, edge contact and two-point
%   contact cannot be found due to low precision for a relatively fast
%   speed.
% 
% 
% Options:
%  - R: the radius of the peg and the hole (they are assumed to be the same size)
%  - d_ecc/h_ecc: parameters for eccentric pegs, currently not used
%  - theta_[x,y,z]: the angle that the peg spins around [x,y,z]-axis
%
% Other settings:
%  - [m,n]_div: used to meshgrid the simulation area
%  - [x,y]: used to set simulation ranges
%
% Author:
%  Ray Lee(lirui2013@ia.ac.cn)
%
% Date:
%  2014-06-01 File created
%  2014-07-17 Add comments; minor structure changed
%  2014-09-14 Create this version on the basis of script_main.m
%  2014-09-23 Add comments

%% Clear the workspace
clear all
close all
clc

%% Set the options
% set simulation parameters
options.R = 3;
options.d_ecc = 1;
options.h_ecc = 3;
options.theta_x = pi/3;
options.theta_y = 0;
options.theta_z = 0;

% Calculate the angle between axis Z_peg and plane XOY(hole) - alpha
% And then the script will calculate ARIE which rotates around Y_hole about alpha
% Finally, the whole figure will be rotated around Z_hole about phi, which
% is equivalant to rotate the peg around X_hole then Y_hole
z = rot([0 1 0], options.theta_y)*rot([1 0 0], options.theta_x)*[0 0 1 1]';
zp = [z(1:2) ; 0];
phi = acos(zp' * [1 0 0]'/norm(zp));
alpha = acos(zp'*z(1:3)/(norm(z(1:3))*norm(zp)));

options.theta_x = 0;
options.theta_y = alpha;

% set simulation precision
m_div = 241;
n_div = 241;
x = linspace(-20,20,m_div);
y = linspace(-20,20,n_div);
[X, Y] = meshgrid(x,y);
d_lowest = zeros(n_div, m_div) + NaN;

% statistics
inters = [0 0 0 0 0];
contacts = [0 0 0 0 0];

% start calculation
for i = 1:length(x)
    for j = 1:length(y)
        param = [x(i) y(j) 5];
%         disp(param);
        
        % 2015
%         [flag, lowest_point, ~] = ARIE_findLowestPoint(param, options, 0);

        [flag, ~, O_pe_H] = ARIE_findLowestPoint(param, options, 0);
        lowest_point = O_pe_H(3);
        
        if flag.outrange == 0 && isnan(lowest_point)
            ARIE_log(param, options, flag);
        else
            inters(flag.inter+1) = inters(flag.inter+1) + 1;
            contacts(flag.contact+2) = contacts(flag.contact+2) + 1;
            clc
            disp('intersection points:   ');
            disp('    0       |     1     |     2     |     3     |     4     ');
            disp(['  ' num2str(inters,'%12.8d')]);
            disp('contact points:   ');
            disp('   edge     |     0     |     1     |     2     |     3     ');
            disp(['  ' num2str(contacts,'%12.8d')]);
        end
        d_lowest(j, i) = lowest_point; 
    end
end

% Rotate back
X1 = X*cos(phi) - Y*sin(phi);
Y1 = X*sin(phi) + Y*cos(phi);

mesh(X1,Y1,d_lowest);

% annotations
xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)');
grid on, axis equal
% axis auto

%% Find the lowest point
z_lowest = min(min(d_lowest));
[row, col] = find(d_lowest == z_lowest);
x_lowest = X1(row,col);
y_lowest = Y1(row,col);
% title(['X = ' num2str(x_lowest(1)) ' Y = ' num2str(y_lowest(1)) ' h = ' num2str(z_lowest)]);
title('Attractive Region in Environment (ARIE) of Peg-Hole System');
fprintf('\n\n');
disp(['Lowest point at: X = '  num2str(x_lowest(1)) ' Y = ' num2str(y_lowest(1))  ' h = ' num2str(z_lowest)]);

%% export the data to .txt file
% ARIE_output(x, y, d_lowest);
