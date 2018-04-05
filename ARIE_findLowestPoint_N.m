function [flag, d_lowest, O_pe_H] = ARIE_findLowestPoint_N(param, options, display)
%ARIE_findLowestPoint_N   Find lowest point of an Attractive Region in 
%                         Environment (Numerical version).
% [INPUT]
%   param:   1x3 vector, the center point of the bottom surface of the peg
%   options: struct, which contains the radius of the peg (and hole), the
%            rotation angles, etc.
%   display: 0 - turn off debug info
%            1 - turn on text debug info
%            2 - turn on figure debug info
%            3 - turn on text & figure debug info
% [OUTPUT]
%   flag:    struct, which stores the status of the intersections and
%            contact points
%   d_lowest:z-coord of the lowest point on the lower surface of the peg
%   O_pe_H:  coord of the center point of the lower surface of the peg 
%
%  Author:
%   Ray Lee(lirui2013@ia.ac.cn)
%
%  LOG:
%   2014-04-11: Create the file.
%   2014-05-20: Fix some small bugs, which solves some display issues in
%               32-bit OS.
%   2014-06-01: Add struct STATISTIC, which is used to count different types
%               of points.
%   2014-07-17: Add output: 'O_pe_H', which corresponds to the center point 
%               of the lower surface of the peg when the peg and hole are
%               touched each other.
%               Fixed a bug on edge_point determination.
%   2014-08-22: Improvements on edge-only detect
%   2015-10-09: Changed to numerical method

% 标记位
flag = [];
flag.outrange = 0;          % 最低点是否在hole范围内，是则为0，不是则为1
flag.inter = 0;             % 投影与底圆交点情况，值为n则有n个交点
flag.contact = 0;           % 轴与孔的接触点情况，值为n则有n个接触点，-1为侧棱接触
% 存放最低点的z坐标值
% d_lowest = NaN;
% 平移量
% x与y为变化量，z为轴相对于孔所在平面抬起高度，设为一个定值即可，比如5
O_x_pe_H = param(1);        % x
O_y_pe_H = param(2);        % y
O_z_pe_H = param(3);        % z

%% 参数
% 轴孔半径
R = options.R;
% 偏心距(沿X_pe轴方向)
% d_ecc = options.d_ecc;
% 偏心距(沿Z_pe方向)
h_ecc = options.h_ecc;

% peg-e坐标轴与hole坐标轴所成夹角
theta_x_pe_H = options.theta_x;
theta_y_pe_H = options.theta_y;
% theta_z_pe_H = options.theta_z;

%% 基准圆，用于各种变换和画图
theta = linspace(0, 2*pi, 721);
rr = R;
xx = rr * cos(theta);
yy = rr * sin(theta);
zz = zeros(1,length(theta));
ww = ones(1,length(theta)); % 齐次坐标

%%
% hole坐标系(基坐标系)------------------------------------------------------
% frame_H = eye(4);
if display == 2 || display == 3
    % 画出hole
    plot3(xx, yy, zz, 'k'); hold on
    % 轴坐标标注
    xlabel('x'); ylabel('y'); zlabel('z');

    grid on; axis equal
    axis([-2*R 2*R -2*R 2*R -2 8]);
end

% peg-e坐标系(偏心轴下部固定)-----------------------------------------------
% peg-e与hole坐标系的关系参数[Xpe Ype Zpe theta_xe theta_ye theta_ze]
axis_o_pe = [1 0 0]';
axis_n_pe = [0 1 0]';
axis_a_pe = [0 0 1]';
O_pe = [0 0 0]';
frame_pe = [axis_o_pe axis_n_pe axis_a_pe O_pe;0 0 0 1];
% peg-e坐标系在hole坐标系下的表示(变换)
T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H);
frame_pe_H = T_pe_H * frame_pe;
% peg-e坐标系的坐标轴、原点在hole坐标系下的表示(变换结果)
% axis_o_pe_H = frame_pe_H(1:3,1);
% axis_n_pe_H = frame_pe_H(1:3,2);
axis_a_pe_H = frame_pe_H(1:3,3);
O_pe_H = frame_pe_H(1:3,4);

%% 画peg-e轴底面圆
% peg-e轴底面圆
points = T_pe_H * [xx;yy;zz;ww];
if display == 2 || display == 3
    plot3(points(1,:),points(2,:),points(3,:),'b');
    hold on
    % peg-e轴底面圆的原点
    plot3(O_pe_H(1),O_pe_H(2),O_pe_H(3),'r.');
    % peg-e轴底面圆在XOY(hole)平面上的投影
    plot3(points(1,:),points(2,:),zz,'b');
end
%% 求最低点

% 20140716-----------------------------------------------------------------
axis_a_pe_H_proj = [axis_a_pe_H(1:2); 0]; % z轴在XoY平面的投影
costt = (axis_a_pe_H_proj' * axis_a_pe_H)/(norm(axis_a_pe_H_proj)*norm(axis_a_pe_H)); %z轴与XoY平面的夹角

vector_R = [0 0 0]'; % 由中心点到最低点的向量
vector_R(3) =  -R * costt; % z坐标由三角形求出
a_R = axis_a_pe_H_proj / norm(axis_a_pe_H_proj) * R * sqrt(1-costt^2);
vector_R(1:2) = a_R(1:2); % x y坐标由投影向量求出

lowest_point = O_pe_H + vector_R;

if display == 2 || display == 3
    % 画轴的底面圆在Z_h方向上的最低点
    plot3([lowest_point(1) lowest_point(1)],...
          [lowest_point(2) lowest_point(2)],...
          [lowest_point(3) 0], 'r');
    hold on
    plot3([lowest_point(1) lowest_point(1)],...
          [lowest_point(2) lowest_point(2)],...
          [lowest_point(3) 0], 'ro');
end

%% 排除最低点投影不在hole内的情况
if norm(lowest_point(1:2)) > R
	if (display == 1 || display == 3), disp('lowest point OUT of range!'), end
    d_lowest = NaN;
    O_pe_H(3) = O_pe_H(3) - lowest_point(3);
    
    flag.outrange = 1;
    return;
end

%% 求圆和椭圆交点

% 2015-09-28 数值方法求圆和椭圆交点
P = create_ellipse_proj(R, theta, T_pe_H); % 求旋转后的peg底面点集
xe = P(1,:);
ye = P(2,:);
f = abs(xe.^2 + ye.^2 - R^2); % 将椭圆上的点带入圆的方程，求极小值点

local_min_idx = localMaximum(-f);
local_min_n = length(local_min_idx);
local_min = theta(local_min_idx);

delta = pi/18; % 10 degree
n_step = 100;
k=1;
while(k < 10) % 迭代10次基本可求出近似解
    delta = delta / 10;
    g=zeros(local_min_n,n_step);
    g_idx=zeros(1,local_min_n);
    for i=1:local_min_n
        pointer = local_min(i);
        
        tt = linspace(pointer-delta, pointer+delta, n_step);
        P = create_ellipse_proj(R, tt, T_pe_H);
        
        xe = P(1,:);
        ye = P(2,:);
        
        g(i,:) = abs(xe.^2 + ye.^2 - R^2);
        t_idx = localMaximum(-g(i,:));
        g_idx(i) = t_idx(1); % 防止t_idx在同一位置求出两个解
        
        local_min(i) = tt(g_idx(i));        
    end
    k = k+1;
end

P = create_ellipse_proj(R, local_min, T_pe_H);
xe = P(1,:);
ye = P(2,:);

f = abs(xe.^2 + ye.^2 - R^2);

% 删去不是交点或切点的点
f_idx = f > 10^-6;
xe(f_idx)=[];
ye(f_idx)=[];

% 首尾求出的可能是同一个点
if (xe(1)-xe(end))^2+(ye(1)-ye(end))^2 < 10^-6
    xe(end)=[];
    ye(end)=[];
end

x = xe;
y = ye;

%% 判断交点情况
% 投影圆和hole的交点有0个、1个、2个、3个、4个5种情况，对应的
% 接触点有以下几种情况：
% 	1. 无接触点                        (0个交点)
% 	2. 有一个接触点                    
% 	2.1  轴的底面圆与孔相切             (1个交点/3个交点)
% 	2.2  轴的底面圆与孔不相切           (2个交点)
%   2.3  轴的侧棱与孔有一个接触点       (2个交点)
% 	3. 有两个接触点                    (2个交点/4个交点)
% 	4. 有三个接触点                    (2个交点/4个交点)

solution_count = length(x);
flag.inter = solution_count;
if solution_count > 4 % 两个圆重合，可以完成insertion
    d_lowest = -h_ecc;
    flag.inter = Inf;
    return;
end


% 求出每个投影交点对应的peg底面上的原象的高度z，通过比较高度判断实际接触点有几个
z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
[~, idx] = sort(z);
z_cmp = abs(z(idx) - z(idx(1)));

switch(solution_count)
    case 0
        d_lowest = 0;
        flag.contact = 0;
        return;
    case 1
        flag.contact = 1;
    case 2
        if norm(z_cmp) > 10^-6 % 两个交点一高一低 - 只有一个接触点
            z = z(idx(1));
            x = x(idx(1));
            y = y(idx(1));
            flag.contact = 1;
        else
            flag.contact = 2;
        end
    case 3
        if norm(z_cmp(1:2)) < 10^-6 % 两个交点一样高 - 有两个接触点
            z = z(idx(1:2));
            x = x(idx(1:2));
            y = y(idx(1:2));
            flag.contact = 2;
        else
            z = z(idx(1));
            x = x(idx(1));
            y = y(idx(1));
            flag.contact = 1;
        end
    case 4
        if norm(z_cmp) < 10^-6 % 四个交点一样高 - 卡在口上装不进去
            d_lowest = 0;
            flag.contact = 4;
            return;
        elseif norm(z_cmp(1:2)) < 10^-6 % 两个交点一样高 - 有两个接触点
            z = z(idx(1:2));
            x = x(idx(1:2));
            y = y(idx(1:2));
            flag.contact = 2;
        else
            z = z(idx(1));
            x = x(idx(1));
            y = y(idx(1));
            flag.contact = 1;
        end
end

%% 求侧棱与孔是否有接触点

% peg底面圆上的点
px = points(1,:);
py = points(2,:);
pz = points(3,:);

% 求peg与三个坐标轴方向的夹角cos(alpha) cos(beta) cos(gamma)
ca = axis_o_pe' * axis_a_pe_H;
cb = axis_n_pe' * axis_a_pe_H;
cc = axis_a_pe' * axis_a_pe_H;
cs = ca^2 + cb^2;

% peg柱面参数方程
% x = x' + t*cos(alpha)
% y = y' + t*cos(beta)
% z = z' + t*cos(gamma)
% t是参数，(x',y',z')是柱面上任意一个圆截面上的点，这里(x',y',z')在peg底面上取点
cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);

% 去除令判别式小于零的点和对应的值
delta_idx = cdelta<0;
px(delta_idx) = NaN;
py(delta_idx) = NaN;
pz(delta_idx) = NaN;
cdelta(delta_idx) = NaN;

% 求出相应的参数t的取值.这时以上参数方程实际上是两圆柱的交线的参数方程
t = [(-(px*ca+py*cb) + sqrt(cdelta)) / cs...
    (-(px*ca+py*cb) - sqrt(cdelta)) / cs];

% 求出交线上的点.当参数t<0时对应的是peg底面以下的交线，实际情况中不存在，故把相
% 应的点和对应的值除去
% ppx = [px px] + t * ca;
% ppy = [py py] + t * cb;
ppz = [pz pz] + t * cc;

t_idx = t<0;
% ppx(t_idx) = NaN;
% ppy(t_idx) = NaN;
ppz(t_idx) = NaN;

[~, idx] = min(ppz);
if idx <= length(t)/2
    delta_flag = 1; % 决定t的解的符号的标志位
else
    delta_flag = 2;
end

thetatheta = [theta theta]; % t的解有两支,对应原向量长度要翻倍
theta_new = thetatheta(idx);

delta = pi/18;
for k=1:20
    delta = delta/2;
    
    theta_aa = theta_new-delta; if theta_aa < 0 , theta_aa = 0; end
    theta_bb = theta_new+delta; if theta_bb > 2*pi , theta_bb = 2*pi; end
    theta_range = linspace(theta_aa, theta_bb, 101);
    P = create_ellipse_proj(R, theta_range, T_pe_H);
    px = P(1,:);
    py = P(2,:);
    pz = P(3,:);
    
    cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);
    
    delta_idx = cdelta<0;
    px(delta_idx) = NaN;
    py(delta_idx) = NaN;
    pz(delta_idx) = NaN;
    cdelta(delta_idx) = NaN;
    
    if delta_flag == 1
        t = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;
    elseif delta_flag == 2
        t = (-(px*ca+py*cb) - sqrt(cdelta)) / cs;
    else
        t = t * NaN;
    end
    
    t_idx = t<0;
    px(t_idx) = NaN;
    py(t_idx) = NaN;
    pz(t_idx) = NaN;
    t(t_idx) = NaN;
    
    ppx = px + t * ca;
    ppy = py + t * cb;
    ppz = pz + t * cc;
    
    [~,min_idx] = min(ppz);
    
    theta_new = theta_range(min_idx);
    edge_min_point = [ppx(min_idx), ppy(min_idx), ppz(min_idx)].';
end

% 比较一下侧棱点和底面点的高度
% edge_x = edge_min_point(1);
% edge_y = edge_min_point(2);
edge_z = edge_min_point(3);

contact_point = [x;y;z];

if abs(edge_z - z(1)) > 10^-6 % 侧面点和底面点高度不同
    if edge_z < z(1) % 侧面接触
        h_down = edge_z;
        contact_point = edge_min_point;
        flag.contact = -1;
    else % 底面接触
        h_down = z(1);
    end
else % 侧面和底面高度相同
    if flag.contact == 2 && abs(t(min_idx)) > 10^-6 % 三点接触,否则底面和侧面是同一个点
        contact_point = [contact_point edge_min_point];
        flag.contact = 3;               
    end
    h_down = z(1);
end

d_lowest = lowest_point(3) - h_down;
O_pe_H(3) = O_pe_H(3) - h_down;

%% Add title to the figure
if display == 2 || display == 3
    plot3(contact_point(1,:), contact_point(2,:), contact_point(3,:), 'r*');
    title(['d_{lowest} = ' num2str(d_lowest)]);
    disp(flag);
end

