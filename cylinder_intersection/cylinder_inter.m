%  Author:
%   Rui Li (raysworld@outlook.com)
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
%   2014-07-18��Distinguish radius of peg and radius of hole
%

%%
clear all
close all
clc
%% ����
% ��ײ���
R = 3;      % ��װ뾶
R_peg = 3.0;
R_hole = 3.0;
d_ecc = 1;  % ƫ�ľ�(��X_pe�᷽��)
h_ecc = 3;  % ƫ�ľ�(��Z_pe����)

% peg-e��hole����ϵ�Ĺ�ϵ����[Xpe Ype Zpe theta_xe theta_ye theta_ze]
theta_x_pe_H = pi/6;       % ���ƿ�(H)����ϵx��ת���ĽǶ�
theta_y_pe_H = 0;         % ���ƿ�(H)����ϵy��ת���ĽǶ�
theta_z_pe_H = 0;                % ���ƿ�(H)����ϵz��ת���ĽǶ�(���趨����Ϊ����Ҫ�ڸ��᷽����ת)
O_x_pe_H = 0;                         % ��ĵ������ĵ��ڿ�����ϵx�᷽��ƫ��׵����ĵ�ľ���
O_y_pe_H = 1;                         % ��ĵ������ĵ��ڿ�����ϵy�᷽��ƫ��׵����ĵ�ľ���
O_z_pe_H = 3;                         % ��ĵ������ĵ��ڿ�����ϵz�᷽��ƫ��׵����ĵ�ľ���

% ����ȫ�������棬�������ݲ����޸�!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%--------------------------------------------------------------------------
%% ��־λ
flag = [];
flag.outrange = 0;          % ��͵��Ƿ���hole��Χ�ڣ�����Ϊ0��������Ϊ1
flag.inter = 0;             % ͶӰ���Բ���������ֵΪn����n������
flag.contact = 0;           % ����׵ĽӴ��������ֵΪn����n���Ӵ��㣬-1Ϊ����Ӵ�
% �����͵��z����ֵ
d_lowest = NaN;

%% ��׼Բ�����ڸ��ֱ任�ͻ�ͼ
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
% hole����ϵ(������ϵ)------------------------------------------------------
frame_H = eye(4);
% ����hole
plot3(xx_h, yy_h, zz_h, 'k');
hold on
% �������ע
xlabel('x'); ylabel('y'); zlabel('z');

grid on
axis equal
axis auto

% peg-e����ϵ(ƫ�����²��̶�)-----------------------------------------------

axis_o_pe = [1 0 0]';
axis_n_pe = [0 1 0]';
axis_a_pe = [0 0 1]';
O_pe = [0 0 0]';
frame_pe = [axis_o_pe axis_n_pe axis_a_pe O_pe;0 0 0 1];
% peg-e����ϵ��hole����ϵ�µı�ʾ(�任)
T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H);
frame_pe_H = T_pe_H * frame_pe;
% peg-e����ϵ�������ᡢԭ����hole����ϵ�µı�ʾ(�任���)
axis_o_pe_H = frame_pe_H(1:3,1);      % peg-e����ϵ��x��
axis_n_pe_H = frame_pe_H(1:3,2);      % peg-e����ϵ��y��
axis_a_pe_H = frame_pe_H(1:3,3);      % peg-e����ϵ��z��
O_pe_H = frame_pe_H(1:3,4);                 % peg-e����ϵ��ԭ��

%% ����peg-e����Բ
% peg-e�����Բ
points = T_pe_H * [xx_p;yy_p;zz_p;ww_p];
plot3(points(1,:),points(2,:),points(3,:),'b');
hold on
% peg-e�����Բ��ԭ��
plot3(O_pe_H(1),O_pe_H(2),O_pe_H(3),'r.');
% peg-e�����Բ��XOY(hole)ƽ���ϵ�ͶӰ
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