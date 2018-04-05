close

%% 轴参数
d = 2;
r = 3;
h = 10;

D = d + r * 2; % 轴心距

%% 基准圆，用于各种变换和画图
t = linspace(0, 2*pi, 360);
rr = r;
xx = rr * cos(t);
yy = rr * sin(t);
zz = zeros(1,length(t));
ww = ones(1,length(t)); % 齐次坐标

% 轴p1,底面圆pb1圆心位于(-(r+d/2),0,0)处
rr_p1 = rr;
xx_p1 = xx - (r+d/2);
yy_p1 = yy;
zz_p1_b = zz;                       % 底面圆z=0
zz_p1_t = h * ones(1,length(t));    % 顶面圆z=h
ww_p1 = ww;

% 轴p2,底面圆pb2圆心位于((r+d/2),0,0)处
rr_p2 = rr;
xx_p2 = xx + (r+d/2);
yy_p2 = yy;
zz_p2_b = zz;                       % 底面圆z=0
zz_p2_t = h * ones(1,length(t));    % 顶面圆z=h
ww_p2 = ww;

plot3(xx_p1,yy_p1,zz_p1_b,'r');hold on
plot3(xx_p1,yy_p1,zz_p1_t,'r');
plot3(xx_p2,yy_p2,zz_p2_b,'b');
plot3(xx_p2,yy_p2,zz_p2_t,'b');

axis equal
grid on
xlabel('x');ylabel('y');zlabel('z');

axis([-15 15 -15 15 -15 15])

%% 旋转
axis_o_pe = [1 0 0]';
axis_n_pe = [0 1 0]';
axis_a_pe = [0 0 1]';
O_pe = [0 0 0]';
frame_pe = [axis_o_pe axis_n_pe axis_a_pe O_pe;0 0 0 1];

% peg-e坐标系在hole坐标系下的表示(变换)
theta_x_pe_H = pi/180*30;          % 轴绕H坐标系x轴转过的角度
theta_y_pe_H = pi/180*45;          % 轴绕H坐标系y轴转过的角度
theta_z_pe_H = pi/180*10;          % 轴绕H坐标系z轴转过的角度

O_x_pe_H = 0;              % 轴的底面中心点在H坐标系x轴方向偏离原点的距离
O_y_pe_H = 0;              % 轴的底面中心点在H坐标系y轴方向偏离原点的距离
O_z_pe_H = 0;              % 轴的底面中心点在H坐标系z轴方向偏离原点的距离
T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H)*rot([0 0 1], theta_z_pe_H);

% 变换后初始双轴各个圆
points_pb1 = T_pe_H * [xx_p1;yy_p1;zz_p1_b;ww_p1];
points_pt1 = T_pe_H * [xx_p1;yy_p1;zz_p1_t;ww_p1];
points_pb2 = T_pe_H * [xx_p2;yy_p2;zz_p2_b;ww_p2];
points_pt2 = T_pe_H * [xx_p2;yy_p2;zz_p2_t;ww_p2];

plot3(points_pb1(1,:),points_pb1(2,:),points_pb1(3,:),'r');hold on
plot3(points_pt1(1,:),points_pt1(2,:),points_pt1(3,:),'r');
plot3(points_pb2(1,:),points_pb2(2,:),points_pb2(3,:),'b');
plot3(points_pt2(1,:),points_pt2(2,:),points_pt2(3,:),'b');

%% 

