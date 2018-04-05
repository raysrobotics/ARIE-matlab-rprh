% 注:此版本输出为底面中心点，不是最低点，其他无区别
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
clear all
close all
clc

%% Set the options
% set simulation parameters
options.R = 3;
options.d_ecc = 1;
options.h_ecc = 3;
options.theta_x = pi/6;
options.theta_y = pi/12;
options.theta_z = 0;

% set simulation precision
m_div = 300;
n_div = 300;
x = linspace(-6,6,m_div);
y = linspace(-6,6,n_div);
[X, Y] = meshgrid(x,y);
d_lowest = zeros(n_div, m_div) + NaN;

% set export settings (file name of the export file)
curr_time = clock;
time_str = [num2str(curr_time(1)) '-' num2str(curr_time(2)) '-' ...
            num2str(curr_time(3)) '-ARIE_OUTPUT_O'];
log_filename = [time_str '.txt'];        
fin = fopen(log_filename, 'a', 'n', 'UTF-8');
format_double = '%4.6f\t';

% start calculation
for i = 1:length(x)
    for j = 1:length(y)
        param = [x(i) y(j) 5];
        [flag, lowest_point, O_point] = ARIE_findLowestPoint(param, options, 0);
                
        if flag.outrange == 0 && isnan(lowest_point)
            continue;
        else
            d_lowest(j, i) = lowest_point;            
            if ~isnan(lowest_point)
                s = [num2str(O_point(1),format_double) '   '...
                          num2str(O_point(2),format_double) '   '...
                          num2str(O_point(3),format_double)];
                 fprintf(fin,'%s\n', s);
            end
        end
    end
end
fclose(fin);

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
