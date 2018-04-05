close

%% �����
d = 2;
r = 3;
h = 10;

D = d + r * 2; % ���ľ�

%% ��׼Բ�����ڸ��ֱ任�ͻ�ͼ
t = linspace(0, 2*pi, 360);
rr = r;
xx = rr * cos(t);
yy = rr * sin(t);
zz = zeros(1,length(t));
ww = ones(1,length(t)); % �������

% ��p1,����Բpb1Բ��λ��(-(r+d/2),0,0)��
rr_p1 = rr;
xx_p1 = xx - (r+d/2);
yy_p1 = yy;
zz_p1_b = zz;                       % ����Բz=0
zz_p1_t = h * ones(1,length(t));    % ����Բz=h
ww_p1 = ww;

% ��p2,����Բpb2Բ��λ��((r+d/2),0,0)��
rr_p2 = rr;
xx_p2 = xx + (r+d/2);
yy_p2 = yy;
zz_p2_b = zz;                       % ����Բz=0
zz_p2_t = h * ones(1,length(t));    % ����Բz=h
ww_p2 = ww;

plot3(xx_p1,yy_p1,zz_p1_b,'r');hold on
plot3(xx_p1,yy_p1,zz_p1_t,'r');
plot3(xx_p2,yy_p2,zz_p2_b,'b');
plot3(xx_p2,yy_p2,zz_p2_t,'b');

axis equal
grid on
xlabel('x');ylabel('y');zlabel('z');

axis([-15 15 -15 15 -15 15])

%% ��ת
axis_o_pe = [1 0 0]';
axis_n_pe = [0 1 0]';
axis_a_pe = [0 0 1]';
O_pe = [0 0 0]';
frame_pe = [axis_o_pe axis_n_pe axis_a_pe O_pe;0 0 0 1];

% peg-e����ϵ��hole����ϵ�µı�ʾ(�任)
theta_x_pe_H = pi/180*30;          % ����H����ϵx��ת���ĽǶ�
theta_y_pe_H = pi/180*45;          % ����H����ϵy��ת���ĽǶ�
theta_z_pe_H = pi/180*10;          % ����H����ϵz��ת���ĽǶ�

O_x_pe_H = 0;              % ��ĵ������ĵ���H����ϵx�᷽��ƫ��ԭ��ľ���
O_y_pe_H = 0;              % ��ĵ������ĵ���H����ϵy�᷽��ƫ��ԭ��ľ���
O_z_pe_H = 0;              % ��ĵ������ĵ���H����ϵz�᷽��ƫ��ԭ��ľ���
T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H)*rot([0 0 1], theta_z_pe_H);

% �任���ʼ˫�����Բ
points_pb1 = T_pe_H * [xx_p1;yy_p1;zz_p1_b;ww_p1];
points_pt1 = T_pe_H * [xx_p1;yy_p1;zz_p1_t;ww_p1];
points_pb2 = T_pe_H * [xx_p2;yy_p2;zz_p2_b;ww_p2];
points_pt2 = T_pe_H * [xx_p2;yy_p2;zz_p2_t;ww_p2];

plot3(points_pb1(1,:),points_pb1(2,:),points_pb1(3,:),'r');hold on
plot3(points_pt1(1,:),points_pt1(2,:),points_pt1(3,:),'r');
plot3(points_pb2(1,:),points_pb2(2,:),points_pb2(3,:),'b');
plot3(points_pt2(1,:),points_pt2(2,:),points_pt2(3,:),'b');

%% 

