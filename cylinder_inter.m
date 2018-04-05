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
%   2014-07-18：Distinguish radius of peg and radius of hole
%

%%
clear all
close all
clc
%% 参数
% 轴孔参数
R = 3;      % 轴孔半径
R_peg = 3.0;
R_hole = 3.0;
d_ecc = 1;  % 偏心距(沿X_pe轴方向)
h_ecc = 3;  % 偏心距(沿Z_pe方向)

% peg-e与hole坐标系的关系参数[Xpe Ype Zpe theta_xe theta_ye theta_ze]
theta_x_pe_H = pi/6;       % 轴绕孔(H)坐标系x轴转过的角度
theta_y_pe_H = 0;         % 轴绕孔(H)坐标系y轴转过的角度
theta_z_pe_H = 0;                % 轴绕孔(H)坐标系z轴转过的角度(不设定，因为不需要在该轴方向旋转)
O_x_pe_H = 0;                         % 轴的底面中心点在空坐标系x轴方向偏离孔的中心点的距离
O_y_pe_H = 1;                         % 轴的底面中心点在空坐标系y轴方向偏离孔的中心点的距离
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
theta = linspace(0, 2*pi, 361);
% rr = R;
% xx = rr * cos(t);
% yy = rr * sin(t);
% zz = zeros(1,length(t));
% ww = ones(1,length(t));

rr_p = R_peg;
xx_p = rr_p * cos(theta);
yy_p = rr_p * sin(theta);
zz_p = zeros(1,length(theta));
ww_p = ones(1,length(theta));

rr_h = R_hole;
xx_h = rr_p * cos(theta);
yy_h = rr_p * sin(theta);
zz_h = zeros(1,length(theta));
ww_h = ones(1,length(theta));


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


%%
px = points(1,:);
py = points(2,:);
pz = points(3,:);

% ca = cos(theta_x_pe_H);
% cb = cos(theta_y_pe_H);
% cc = cos(theta_z_pe_H);

ca = axis_o_pe' * axis_a_pe_H;
cb = axis_n_pe' * axis_a_pe_H;
cc = axis_a_pe' * axis_a_pe_H;
cs = ca^2 + cb^2;

cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);

idx = cdelta<0;
px(idx) = NaN;
py(idx) = NaN;
pz(idx) = NaN;
cdelta(idx) = NaN;

t1 = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;
t2 = (-(px*ca+py*cb) - sqrt(cdelta)) / cs;
% t=linspace(0,10,292);
%
% ppx = px + t * ca;
% ppy = py + t * cb;
% ppz = pz + t * cc;
%
% plot3(ppx, ppy, ppz,'r.'); hold on

% ppx1 = px + t1 * ca;
% ppy1 = py + t1 * cb;
% ppz1 = pz + t1 * cc;
%
% plot3(ppx1, ppy1, ppz1,'r.'); hold on
%
% ppx2 = px + t2 * ca;
% ppy2 = py + t2 * cb;
% ppz2 = pz + t2 * cc;
%
% plot3(ppx2, ppy2, ppz2,'b.');

tt = [t1 t2];
ppx = [px px] + tt * ca;
ppy = [py py] + tt * cb;
ppz = [pz pz] + tt * cc;

tidx = tt<0;
ppx(tidx) = NaN;
ppy(tidx) = NaN;
ppz(tidx) = NaN;

plot3(ppx, ppy, ppz,'r.'); hold on
%%
[~,idx] = min(ppz);
plot3(ppx(idx), ppy(idx), ppz(idx), 'k*');
tt(idx)
% %%
k=1;
delta = pi/2;
tt_new = [theta theta];
t_new = tt_new(idx);
while(k < 10)
    
    delta = delta/10;
    ttt = linspace(t_new-delta, t_new+delta, 100);
    P = create_ellipse_proj(R, ttt, T_pe_H);
    px = P(1,:);
    py = P(2,:);
    pz = P(3,:);
%     plot3(px,py,pz,'go');
    
    cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);
    
    idx = cdelta<0;
    px(idx) = NaN;
    py(idx) = NaN;
    pz(idx) = NaN;
    cdelta(idx) = NaN;
    
    t1 = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;
    t2 = (-(px*ca+py*cb) - sqrt(cdelta)) / cs;
    
    tt = [t1 t2];
    ppx = [px px] + tt * ca;
    ppy = [py py] + tt * cb;
    ppz = [pz pz] + tt * cc;
    
    tidx = tt<0;
    ppx(tidx) = NaN;
    ppy(tidx) = NaN;
    ppz(tidx) = NaN;
    
    [~,idx] = min(ppz);
    plot3(ppx(idx), ppy(idx), ppz(idx), 'k*');
    tt(idx)
    
    tt_new = [ttt ttt];
    t_new = tt_new(idx);
    
    k=k+1
end