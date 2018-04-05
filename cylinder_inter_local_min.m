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
theta_x_pe_H = pi/4;       % 轴绕孔(H)坐标系x轴转过的角度
theta_y_pe_H = pi/6;         % 轴绕孔(H)坐标系y轴转过的角度
theta_z_pe_H = 0;                % 轴绕孔(H)坐标系z轴转过的角度(不设定，因为不需要在该轴方向旋转)
O_x_pe_H = 0.6;                         % 轴的底面中心点在空坐标系x轴方向偏离孔的中心点的距离
O_y_pe_H = 2.6;                         % 轴的底面中心点在空坐标系y轴方向偏离孔的中心点的距离
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
theta = linspace(0, 2*pi, 3601);
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
plot3(points(1,:),points(2,:),points(3,:),'k');
hold on
% peg-e轴底面圆的原点
plot3(O_pe_H(1),O_pe_H(2),O_pe_H(3),'r.');
% peg-e轴底面圆在XOY(hole)平面上的投影
plot3(points(1,:),points(2,:),zz_p,'k');


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

delta_idx = cdelta<0;
px(delta_idx) = NaN;
py(delta_idx) = NaN;
pz(delta_idx) = NaN;
cdelta(delta_idx) = NaN;

ppx=zeros(2, length(px));
ppy=zeros(2, length(px));
ppz=zeros(2, length(px));

t = [(-(px*ca+py*cb) + sqrt(cdelta)) / cs;...
    (-(px*ca+py*cb) - sqrt(cdelta)) / cs];

figure(4);
subplot(4,1,1);plot(theta,t(1,:));title('t(\Delta>0)');
subplot(4,1,2);plot(theta,t(2,:));title('t(\Delta<0)');
figure(1);

for i=1:2
    ppx(i,:) = px + t(i,:) * ca;
    ppy(i,:) = py + t(i,:) * cb;
    ppz(i,:) = pz + t(i,:) * cc;
    
    t_idx = t(i,:)<0;
    ppx(i,t_idx) = NaN;
    ppy(i,t_idx) = NaN;
    ppz(i,t_idx) = NaN;
    
    plot3(ppx(i,:),ppy(i,:),ppz(i,:),'Marker','.','Color',[rand(1) rand(1) rand(1)]);
    
    figure(4);
    subplot(4,1,i);hold on;plot(theta(t_idx),t(i,t_idx),'r.');
    figure(1);
end



% for i=1:length(ppx(1,:))
%     plot3(ppx(1,i),ppy(1,i),ppz(1,i),'Marker','.','Color','r');
% %     pause(0.05)
% end
% for i=1:length(ppx(1,:))
%     plot3(ppx(2,i),ppy(2,i),ppz(2,i),'Marker','.','Color','b');
% %     pause(0.05)
% end

real_ppz = ppz(1,:);
real_idx = isnan(ppz(1,:));
real_ppz(real_idx) = [];
idx1 = localMaximum(-real_ppz);
for i=1:length(idx1)
    temp = find(ppz(1,:) == real_ppz(idx1(i)));
    idx1(i) = temp(1);
end

real_ppz = ppz(2,:);
real_idx = isnan(ppz(2,:));
real_ppz(real_idx) = [];
idx2 = localMaximum(-real_ppz);
for i=1:length(idx2)
    temp = find(ppz(2,:) == real_ppz(idx2(i)));
    idx2(i) = temp;
end
% idx1 = localMaximum(-ppz(1,:));
% idx2 = localMaximum(-ppz(2,:));
% % 标记初始极值点
% plot3(ppx(1,idx1), ppy(1,idx1), ppz(1,idx1), 'rd');
% plot3(ppx(2,idx2), ppy(2,idx2), ppz(2,idx2), 'bd');

figure(4);
subplot(4,1,3);plot(theta,ppz(1,:));title('z(\Delta>0)');
subplot(4,1,4);plot(theta,ppz(2,:));title('z(\Delta<0)');

subplot(4,1,3);hold on;plot(theta(idx1),ppz(1,idx1), 'g*');
subplot(4,1,4);hold on;plot(theta(idx2),ppz(2,idx2), 'g*');
figure(1);


% for 1
k=1;
delta1 = pi;
delta2 = pi;

theta_new1 = theta(idx1);
tt_new1 = zeros(1, length(idx1));
local_min_point1 = zeros(3, length(idx1));
ttt1 = zeros(length(idx1), 101);

theta_new2 = theta(idx2);
tt_new2 = zeros(1, length(idx2));
local_min_point2 = zeros(3, length(idx2));
ttt2 = zeros(length(idx2), 101);
while(k < 10)
    delta1 = delta1/2;
    delta2 = delta2/2;
    
    for i=1:length(idx1)
        aa = theta_new1(i)-delta1; if aa < 0 , aa = 0; end
        bb = theta_new1(i)+delta1; if bb > 2*pi , bb = 2*pi; end
        ttt1(i,:) = linspace(aa, bb, 101);
        P = create_ellipse_proj(R, ttt1(i,:), T_pe_H);
        px = P(1,:);
        py = P(2,:);
        pz = P(3,:);
        
        cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);
        
        delta_idx = cdelta<0;
        px(delta_idx) = NaN;
        py(delta_idx) = NaN;
        pz(delta_idx) = NaN;
        cdelta(delta_idx) = NaN;
        
        t1 = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;
        
        t_idx = t1<0;
        px(t_idx) = NaN;
        py(t_idx) = NaN;
        pz(t_idx) = NaN;
        t1(t_idx) = NaN;
        
        ppx1 = px + t1 * ca;
        ppy1 = py + t1 * cb;
        ppz1 = pz + t1 * cc;
        
        [~,min_idx] = min(ppz1);
        %         plot3(ppx1(min_idx), ppy1(min_idx), ppz1(min_idx), 'gd');
        
        theta_new1(i) = ttt1(i, min_idx);
        tt_new1(i) = t1(min_idx);
        local_min_point1(:,i) = [ppx1(min_idx), ppy1(min_idx), ppz1(min_idx)].';
        
        figure(4);
        subplot(4,1,1);hold on;plot(theta(min_idx),t(1,min_idx),'g*');
        figure(1);
    end
    %     plot3(local_min_point1(1,:), local_min_point1(2,:), local_min_point1(3,:), 'rd');
    
    for i=1:length(idx2)
        aa = theta_new2(i)-delta2; if aa < 0 , aa = 0; end
        bb = theta_new2(i)+delta2; if bb > 2*pi , bb = 2*pi; end
        ttt2(i,:) = linspace(aa, bb, 101);
        P = create_ellipse_proj(R, ttt2(i,:), T_pe_H);
        px = P(1,:);
        py = P(2,:);
        pz = P(3,:);
        
        cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);
        
        delta_idx = cdelta<0;
        px(delta_idx) = NaN;
        py(delta_idx) = NaN;
        pz(delta_idx) = NaN;
        cdelta(delta_idx) = NaN;
        
        t2 = (-(px*ca+py*cb) - sqrt(cdelta)) / cs;
        
        t_idx = t2<0;
        px(t_idx) = NaN;
        py(t_idx) = NaN;
        pz(t_idx) = NaN;
        t2(t_idx) = NaN;
        
        ppx2 = px + t2 * ca;
        ppy2 = py + t2 * cb;
        ppz2 = pz + t2 * cc;
        
        %         plot3(ppx2, ppy2, ppz2, 'Marker','s','Color',[rand(1) rand(1) rand(1)]); hold on
        
        [~,min_idx] = min(ppz2);
        
        %         plot3(ppx2(min_idx), ppy2(min_idx), ppz2(min_idx), 'gd');
        %         pause;
        
        theta_new2(i) = ttt2(i, min_idx);
        tt_new2(i) = t2(min_idx);
        local_min_point2(:,i) = [ppx2(min_idx), ppy2(min_idx), ppz2(min_idx)].';
   
        figure(4);
        subplot(4,1,2);hold on;plot(theta(min_idx),t(2,min_idx),'g*');
        figure(1);
    end
    
    k=k+1;
    %     local_min_point1
    %     local_min_point2
    %     plot3(local_min_point1(1,:), local_min_point1(2,:), local_min_point1(3,:), 'gd');
    %     plot3(local_min_point2(1,:), local_min_point2(2,:), local_min_point2(3,:), 'gs');
    %     pause(0.5);
    local_min_point1
    local_min_point2
end
% local_min_point1
% local_min_point2
plot3(local_min_point1(1,:), local_min_point1(2,:), local_min_point1(3,:), 'r*');
plot3(local_min_point2(1,:), local_min_point2(2,:), local_min_point2(3,:), 'b*');

figure(4)
subplot(4,1,4);axis([0 7 0 10]);