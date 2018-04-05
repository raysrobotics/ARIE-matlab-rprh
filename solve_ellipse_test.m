%%
clear
% close all
clc

%%
theta_x_pe_H = 0.39794;       % 轴绕孔(H)坐标系x轴转过的角度
theta_y_pe_H = 0.10472;         % 轴绕孔(H)坐标系y轴转过的角度
theta_z_pe_H = 0;                % 轴绕孔(H)坐标系z轴转过的角度(不设定，因为不需要在该轴方向旋转)
O_x_pe_H = 0;                         % 轴的底面中心点在空坐标系x轴方向偏离孔的中心点的距离
O_y_pe_H = 0;                         % 轴的底面中心点在空坐标系y轴方向偏离孔的中心点的距离
O_z_pe_H = 5;                         % 轴的底面中心点在空坐标系z轴方向偏离孔的中心点的距离

%%
n_div = 300;
tt = linspace(0,2*pi,n_div);
rr = 3;
xx = rr*cos(tt);
yy = rr*sin(tt);
zz = 0*tt;
ww = 0*tt+1;

T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H);
P = create_ellipse_proj(rr, tt, T_pe_H);

xe = P(1,:);
ye = P(2,:);

plot(xe,ye); hold on
grid on
axis equal

xr = rr*cos(tt);
yr = rr*sin(tt);
plot(xr,yr);
%%
f = abs(xe.^2 + ye.^2 - rr^2);
figure;plot(tt, f)
grid on

local_min_idx = localMaximum(-f);
local_min_n = length(local_min_idx);
local_min = tt(local_min_idx);

delta = pi/10; % 5 degree
n_step = 100;
k=1;
while(k < 10)
    delta = delta / 10;
    f=zeros(local_min_n,n_step);
    f_idx=zeros(1,local_min_n);
    for i=1:local_min_n
        pointer = local_min(i);
        
        ttt = linspace(pointer-delta, pointer+delta, n_step);
        P = create_ellipse_proj(rr, ttt, T_pe_H);
        
        xe = P(1,:);
        ye = P(2,:);
        
        f(i,:) = abs(xe.^2 + ye.^2 - rr^2);
        t_idx = localMaximum(-f(i,:));
        f_idx(i) = t_idx(1); % 防止t_idx在同一位置求出两个解
        
        local_min(i) = ttt(f_idx(i));
        %         figure;plot(ttt, f(i,:));hold on;plot(ttt(f_idx(i)), f(i,f_idx(i)), '*');
    end
    P = create_ellipse_proj(rr, local_min, T_pe_H);
    xe = P(1,:);
    ye = P(2,:);
    
    f = abs(xe.^2 + ye.^2 - rr^2);
    
    k = k+1;
end

% 删去不是交点或切点的点
f_idx = f > 10^-6;
xe(f_idx)=[];
ye(f_idx)=[];
f(f_idx) = [];
    disp(xe); disp(ye); disp(f);
