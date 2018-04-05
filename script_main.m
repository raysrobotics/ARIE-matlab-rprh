% This .m script is used for visualizing the ARIE 
% in the assemby task of round-peg and round-hole
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

%% Clear the workspace
clear
close all
clc

%% Set the options
% set simulation parameters
options.R = 3;
options.d_ecc = 1;
options.h_ecc = 3;
options.theta_x = 0;
options.theta_y = pi/4;
options.theta_z = 0;

% set simulation precision
m_div = 241;
n_div = 241;
x = linspace(-8,8,m_div);
y = linspace(-4,4,n_div);
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
        [flag, lowest_point, ~] = ARIE_findLowestPoint(param, options, 0);
        
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
mesh(X,Y,d_lowest);

% annotations
xlabel('x'); ylabel('y'); zlabel('z');
grid on, axis equal
% axis auto

%% Find the lowest point
z_lowest = min(min(d_lowest));
[row, col] = find(d_lowest == z_lowest);
x_lowest = X(row,col);
y_lowest = Y(row,col);
title(['X = ' num2str(x_lowest(1)) ' Y = ' num2str(y_lowest(1)) ' h = ' num2str(z_lowest)]);
fprintf('\n\n');
disp(['Lowest point at: X = '  num2str(x_lowest(1)) ' Y = ' num2str(y_lowest(1))  ' h = ' num2str(z_lowest)]);

%% export the data to .txt file
% ARIE_output(x, y, d_lowest);
