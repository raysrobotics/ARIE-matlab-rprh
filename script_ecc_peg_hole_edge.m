%  Author:
%   Ray Lee(lirui2013@ia.ac.cn)
%
%  LOG:
%   2014-04-11: Create the file.
%   2014-05-20: Fix some small bugs, which solves some display issues in
%               32-bit OS.
%   2014-07-16: Add some comments
%   2014-07-17: Add output: 'O_pe_H', which corresponds to the center point 
%               of the lower surface of the peg when the peg and hole are
%               touched each other.
%               Fixed a bug on edge_point determination.
%   2014-07-18: Distinguish radius of peg and radius of hole
%   2014-08-22: Improvements on edge-only detect
%

%% 
clear all
close all
clc
%% 参数
% 轴孔参数
R = 30;      % 轴孔半径
R_peg = 30.0;
R_hole = 30.2;
d_ecc = 1;  % 偏心距(沿X_pe轴方向)
h_ecc = 3;  % 偏心距(沿Z_pe方向)

% peg-e与hole坐标系的关系参数[Xpe Ype Zpe theta_xe theta_ye theta_ze]
theta_x_pe_H = pi/3;       % 轴绕孔(H)坐标系x轴转过的角度
theta_y_pe_H = 0;         % 轴绕孔(H)坐标系y轴转过的角度
theta_z_pe_H = 0;                % 轴绕孔(H)坐标系z轴转过的角度(不设定，因为不需要在该轴方向旋转)
O_x_pe_H = -0.4667;                         % 轴的底面中心点在空坐标系x轴方向偏离孔的中心点的距离
O_y_pe_H = 0.125;                         % 轴的底面中心点在空坐标系y轴方向偏离孔的中心点的距离
O_z_pe_H = 3;                         % 轴的底面中心点在空坐标系z轴方向偏离孔的中心点的距离

% 参数全部在上面，以下内容不需修改!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%--------------------------------------------------------------------------
%% 标志位
flag = [];
flag.outrange = 0;          % 最低点是否在hole范围内，是则为0，不是则为1
flag.inter = 0;             % 投影与底圆交点情况，值为n则有n个交点
flag.contact = 0;           % 轴与孔的接触点情况，值为n则有n个接触点，-1为侧棱接触
% 存放最低点的z坐标值
d_lowest = NaN;

%% 基准圆，用于各种变换和画图
t = linspace(0, 2*pi, 360);
% rr = R;
% xx = rr * cos(t);
% yy = rr * sin(t);
% zz = zeros(1,length(t));
% ww = ones(1,length(t));

rr_p = R_peg;
xx_p = rr_p * cos(t);
yy_p = rr_p * sin(t);
zz_p = zeros(1,length(t));
ww_p = ones(1,length(t));

rr_h = R_hole;
xx_h = rr_p * cos(t);
yy_h = rr_p * sin(t);
zz_h = zeros(1,length(t));
ww_h = ones(1,length(t));


%%
% hole坐标系(基坐标系)------------------------------------------------------
frame_H = eye(4);
% 画出hole
plot3(xx_h, yy_h, zz_h, 'k');
hold on
% 轴坐标标注
xlabel('x'); ylabel('y'); zlabel('z');

grid on
axis equal
axis auto

% peg-e坐标系(偏心轴下部固定)-----------------------------------------------

axis_o_pe = [1 0 0]';
axis_n_pe = [0 1 0]';
axis_a_pe = [0 0 1]';
O_pe = [0 0 0]';
frame_pe = [axis_o_pe axis_n_pe axis_a_pe O_pe;0 0 0 1];
% peg-e坐标系在hole坐标系下的表示(变换)
T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H);
frame_pe_H = T_pe_H * frame_pe;
% peg-e坐标系的坐标轴、原点在hole坐标系下的表示(变换结果)
axis_o_pe_H = frame_pe_H(1:3,1);      % peg-e坐标系的x轴
axis_n_pe_H = frame_pe_H(1:3,2);      % peg-e坐标系的y轴
axis_a_pe_H = frame_pe_H(1:3,3);      % peg-e坐标系的z轴
O_pe_H = frame_pe_H(1:3,4);                 % peg-e坐标系的原点

%% 画出peg-e底面圆
% peg-e轴底面圆
points = T_pe_H * [xx_p;yy_p;zz_p;ww_p];
plot3(points(1,:),points(2,:),points(3,:),'b');
hold on
% peg-e轴底面圆的原点
plot3(O_pe_H(1),O_pe_H(2),O_pe_H(3),'r.');
% peg-e轴底面圆在XOY(hole)平面上的投影
plot3(points(1,:),points(2,:),zz_p,'b');

%% 求最低点

% legacy code-----------------------------------------------------------------
% % 求轴的底面圆在Z_h方向上的最低点
% [~, j]=find(points == min(min(points(3,:))));

% 20140716-----------------------------------------------------------------
axis_a_pe_H_proj = [axis_a_pe_H(1:2); 0]; % z轴在XoY平面的投影
costt = (axis_a_pe_H_proj' * axis_a_pe_H)/(norm(axis_a_pe_H_proj)*norm(axis_a_pe_H)); %z轴与XoY平面的夹角

vector_R = [0 0 0]'; % 由中心点到最低点的向量
vector_R(3) =  -R_hole * costt; % z坐标由三角形求出
a_R = axis_a_pe_H_proj / norm(axis_a_pe_H_proj) * R_hole * sqrt(1-costt^2);
vector_R(1:2) = a_R(1:2); % x y坐标由投影向量求出

lowest_point = O_pe_H + vector_R;

% 画轴的底面圆在Z_h方向上的最低点
plot3([lowest_point(1) lowest_point(1)],...
      [lowest_point(2) lowest_point(2)],...
      [lowest_point(3) 0], 'r');
hold on
plot3([lowest_point(1) lowest_point(1)],...
      [lowest_point(2) lowest_point(2)],...
      [lowest_point(3) 0], 'ro');
%% 排除最低点投影不在hole内的情况
if norm(lowest_point(1:2)) > R_hole
    disp('lowest point OUT of range!');
	d_lowest = NaN;
	flag.outrange = 1;
end

%%
% 求底面圆和投影圆的交点
% syms x y
% Px = sin(theta_x_pe_H)*tan(theta_x_pe_H) + cos(theta_x_pe_H);
% Py = sin(theta_y_pe_H)*tan(theta_y_pe_H) + cos(theta_y_pe_H);
% Q = tan(theta_x_pe_H)*tan(theta_y_pe_H);

% f1 = x^2 + y^2 - R^2;
% f2 = (x - O_x_pe_H)^2*Py^2 - 2*(x - O_x_pe_H)*(y - O_y_pe_H)*Q + (y - O_y_pe_H)^2*(Q^2+Px^2) - R^2;
% 
% S = solve(f1,f2,'Real',true);
% 
% A = (S.x == real(S.x));
% x = eval(S.x(A));
% B = (S.y == real(S.y));
% y = eval(S.y(B));

% 20140509-----------------------------------------------------------------
% f1 = x^2 + y^2 - R^2;
A1 = 1; B1 = 1; C1 = 0; D1 = 0; E1 = 0; F1 = -R_hole^2;
% f2 = (x - O_x_pe_H)^2*Py^2 - 2*(x - O_x_pe_H)*(y - O_y_pe_H)*Q + (y - O_y_pe_H)^2*(Q^2+Px^2) - R^2;
%A2 = Py^2;
%B2 = Px^2 + Q^2;
%C2 = - Q;
%D2 = Q*O_y_pe_H - Py^2*O_x_pe_H;
%E2 = Q*O_x_pe_H - (Q^2 + Px^2)*O_y_pe_H;
%F2 = Py^2*O_x_pe_H^2 + Px^2*O_y_pe_H^2 + Q^2*O_y_pe_H^2 - 2*Q*O_x_pe_H*O_y_pe_H - R^2;

% 直接使用方程求解存在误差，换为从椭圆上取点，然后求系数的方法
sel_points = [1 73 145 217 289];
coe_xx = points(1,sel_points)';
coe_yy = points(2,sel_points)';

coe_C = [coe_xx.^2, coe_yy.^2, 2*coe_xx.*coe_yy, 2*coe_xx, 2*coe_yy];
sol_C = coe_C\ones(size(coe_xx)); % 通过椭圆上的五个点求A、B、C、D、E五个系数，F=1

A2 = sol_C(1);
B2 = sol_C(2);
C2 = sol_C(3);
D2 = sol_C(4);
E2 = sol_C(5);
F2 = -1;
f1 = [A1 C1 D1; C1 B1 E1; D1 E1 F1];
f2 = [A2 C2 D2; C2 B2 E2; D2 E2 F2];
S = intersectConics(f1, f2);
% plot(S(1,:) ./ S(3,:) , S(2,:) ./ S(3,:), 'ro');
[m_s, ~] = size(S);
if m_s == 3
    x = S(1,:) ./ S(3,:);
    y = S(2,:) ./ S(3,:);
else
    x=[];
    y=[];
end

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
% 先将四个交点的情况处理为两个交点
if solution_count == 4 % 比较四个点的高度，留下高度低的两个
    disp('four intersection points!');
	flag.inter = 4;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    [~, idx] = sort(z);
    z = z(idx(1:2));
    x = x(idx(1:2));
    y = y(idx(1:2));
end
if solution_count == 2 % 比较两个点是否是同一个点(2014-07-18)
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    if norm([x(1)-x(2) y(1)-y(2) z(1)-z(2)]) < 10^-6
        x = x(1);
        y = y(1);
        z = z(1);
    end
end

solution_count = length(x);
% 0个交点
if solution_count == 0
    disp('NO intersection points!');
    flag.inter = 0;
elseif solution_count == 1
    disp('ONE intersection point!');
	flag.inter = 1;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    h_down = z;
    d_lowest = lowest_point(3) - h_down;
    
% 20140718-----------------------------------------------------------------
    disp('ONE contact point!');
    flag.contact = 1;
% 20140718-END-------------------------------------------------------------
elseif solution_count == 3
    disp('THREE intersection points!');
	flag.inter = 3;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    [~, idx] = sort(z);
    z = z(idx(1));
    x = x(idx(1));
    y = y(idx(1));
    
    h_down = z;
    d_lowest = lowest_point(3) - h_down;
elseif solution_count == 2
    if flag.inter ~= 4
        disp('TWO intersection points!');
        flag.inter = 2;
    end
    
    % 计算交点位置到hole平面的距离
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    plot3([x(1) x(1)],[y(1) y(1)],[0 z(1)], 'g-');
    plot3([x(2) x(2)],[y(2) y(2)],[0 z(2)], 'g-');
    
    % 保存交点坐标到新的变量中
    inter_x = x;
    inter_y = y;
    inter_z = z;
    
    if norm(z(1)-z(2)) < 10^-4 % 两接触点或三接触点
        h_down = z(1);
        % 以下判断是否存在第三个接触点
        % 先求一个垂直于XoY平面的平面,该平面过Ope点和最低点
        % 平面方程为Ax+By=1
        if abs(O_pe_H(1)-lowest_point(1))<0.02
            coe_planar = [1,0,O_pe_H(1)];
        elseif abs(O_pe_H(2)-lowest_point(2))<0.02
            coe_planar = [0,1,O_pe_H(2)];
        else
            coe_planar = [O_pe_H(1) O_pe_H(2);lowest_point(1) lowest_point(2)]\[1;1];
            coe_planar(3) = 1;
        end
        
        cy = cos(theta_y_pe_H);
        sy = sin(theta_y_pe_H);
        cx = cos(theta_x_pe_H);
        sx = sin(theta_x_pe_H);
        syms x y
            f11  = (x - O_x_pe_H)*cy-(h_down-O_z_pe_H)*sy;
            f121 = (x - O_x_pe_H)*sx*sy;
            f122 = (y - O_y_pe_H)*cx;
            f123 = (h_down-O_z_pe_H)*sx*cy;
        f1 = f11^2 + (f121 + f122 + f123)^2 - R_peg^2;        
        f2 = coe_planar(1)*x + coe_planar(2)*y - coe_planar(3)*1;
        
        S = solve(f2, f1, 'Real', true);
        
        if isempty(S)
            disp('Equation solve failed!');
            d_lowest = NaN;
            return;
        end
        
        A = (S.x == real(S.x));
        x = eval(S.x(A));
        B = (S.y == real(S.y));
        y = eval(S.y(B));
        
%         plot3(x(1), y(1), h_down, 'go');
%         plot3(x(2), y(2), h_down, 'go');
        
        % 判断哪个解是棱上的点(Ope、最低点、棱上点三点连线成直角)
        v1 = O_pe_H - lowest_point;
        v2 = [x(1);y(1);h_down] - lowest_point;
        if abs(v1.' * v2) < 10^-6
            edge_point = [x(1) y(1) h_down].';
        else
            edge_point = [x(2) y(2) h_down].';
        end
        plot3(edge_point(1), edge_point(2), edge_point(3), 'mo');hold on
        plot3([edge_point(1) lowest_point(1)],...
              [edge_point(2) lowest_point(2)],...
              [edge_point(3) lowest_point(3)], 'm');
        
        % 判断棱上一点的投影与hole的位置关系
        if abs(edge_point(1)^2 + edge_point(2)^2 - R_hole^2) <= 10^-4 % 落在hole上，三点接触
            disp('three-point contact!');
            flag.contact = 3;
            d_lowest = lowest_point(3) - h_down;
        elseif edge_point(1)^2 + edge_point(2)^2 - R_hole^2 < 0 % 落在hole内，两点接触
            disp('two-point contact!');
            flag.contact = 2;
            d_lowest = lowest_point(3) - h_down;
        elseif edge_point(1)^2 + edge_point(2)^2 - R_hole^2 > 0 % 落在hole外，侧棱点接触
            disp('edge-only contact!');
            flag.contact = -1;
            
            d1 = norm(lowest_point(1:2) - edge_point(1:2));
            mid_point = [sum(inter_x)/2  sum(inter_y)/2  sum(inter_z)/2].';
            
            d2_cos = ([0;0] - mid_point(1:2))' * (edge_point(1:2) - mid_point(1:2));
%             d2_temp = norm(mid_point(1:2) - edge_point(1:2)) - norm(mid_point(1:2)) - R_hole;
            if d2_cos >= 0 % mid_point在O_hole和edge_point_proj的同侧
                d2 = norm(mid_point(1:2) - edge_point(1:2)) - norm(mid_point(1:2)) - R_hole;
            else % mid_point在O_hole和edge_point_proj的中间
                d2 = norm(mid_point(1:2) - edge_point(1:2)) + norm(mid_point(1:2)) - R_hole;
            end
                        
            d3 = abs(h_down - lowest_point(3));
            d4 = d2*d3/d1;
            
            d_lowest = lowest_point(3) - (h_down - d4);
        end
    else % 两点高度相差较大，说明只有一个接触点，此处再判断一下是不是侧棱接触(2014-08-22)        
        px = points(1,:);
        py = points(2,:);
        pz = points(3,:);

        ca = axis_o_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cb = axis_n_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cc = axis_a_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cs = ca^2 + cb^2;

        cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);

        idx = cdelta>=0;
        
        px = px(idx);
        py = py(idx);
        pz = pz(idx);
        cdelta = cdelta(idx);

        t1 = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;

        ppx = px + t1 * ca;
        ppy = py + t1 * cb;
        ppz = pz + t1 * cc;

        plot3(ppx, ppy, ppz,'r.');
       
        [~, idx] = min(ppz);
        plot3(ppx(idx), ppy(idx), ppz(idx), 'ko');

        if abs(O_pe_H(1)-lowest_point(1))<10^-6%0.02  % 平面接近与X轴或y轴重合时，方程病态，不能求逆，故直接计算
            coe_planar = [1,0,O_pe_H(1)];
        elseif abs(O_pe_H(2)-lowest_point(2))<10^-6%0.02
            coe_planar = [0,1,O_pe_H(2)];
        else
            coe_planar = [O_pe_H(1) O_pe_H(2);lowest_point(1) lowest_point(2)]\[1;1];
            coe_planar(3) = 1;
        end
        
        cd = abs(coe_planar(1)*ppx + coe_planar(2)*ppy - coe_planar(3)*1) < 10^-4;
        [~, idx] = min(ppz(cd));
        
        if isempty(idx) || min(z) < ppz(idx) % 一点接触            
            disp('one-point contact!');
            flag.contact = 1;

            h_down = min(z);
        else
           disp('edge-only contact!');
            flag.contact = -1;
            
            h_down = ppz(idx);            
        end
        d_lowest = lowest_point(3) - h_down;
    end
end

%% 异常点检测
if d_lowest > 5 || d_lowest < -5
    d_lowest = NaN;
    fprintf(2,'Unexpected result!\n');
end

%% 最低点数值
title(['d_{lowest} = ' num2str(d_lowest)]);
disp(flag);





