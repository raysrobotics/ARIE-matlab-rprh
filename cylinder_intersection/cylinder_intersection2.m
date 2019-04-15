%  ����Բ�����ߵ���͵㣬���жϸõ��Ƿ��ڵ�����
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
theta_x_pe_H = pi/4;       % ���ƿ�(H)����ϵx��ת���ĽǶ�
theta_y_pe_H = pi/4;         % ���ƿ�(H)����ϵy��ת���ĽǶ�
theta_z_pe_H = 0;                % ���ƿ�(H)����ϵz��ת���ĽǶ�(���趨����Ϊ����Ҫ�ڸ��᷽����ת)
O_x_pe_H = -0.6;                         % ��ĵ������ĵ��ڿ�����ϵx�᷽��ƫ��׵����ĵ�ľ���
O_y_pe_H = 1.6;                         % ��ĵ������ĵ��ڿ�����ϵy�᷽��ƫ��׵����ĵ�ľ���
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
theta = linspace(0, 2*pi, 721);

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
plot3(points(1,:),points(2,:),points(3,:),'k');
hold on
% peg-e�����Բ��ԭ��
plot3(O_pe_H(1),O_pe_H(2),O_pe_H(3),'r.');
% peg-e�����Բ��XOY(hole)ƽ���ϵ�ͶӰ
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

t = [(-(px*ca+py*cb) + sqrt(cdelta)) / cs...
    (-(px*ca+py*cb) - sqrt(cdelta)) / cs];

ppx = [px px] + t * ca;
ppy = [py py] + t * cb;
ppz = [pz pz] + t * cc;

t_idx = t<0;
ppx(t_idx) = NaN;
ppy(t_idx) = NaN;
ppz(t_idx) = NaN;

plot3(ppx,ppy,ppz,'Marker','.','Color',[rand(1) rand(1) rand(1)]);

[~, idx] = min(ppz);
if idx <= length(t)/2
    delta_flag = 1; % ����t�Ľ�ķ��ŵı�־λ
else
    delta_flag = 2;
end

thetatheta = [theta theta];
theta_new = thetatheta(idx);
tt_new = 0;
ttt = zeros(1, 101);

k=1;
delta = pi/10;
while(k < 20)
    delta = delta/2;
    
    aa = theta_new-delta; if aa < 0 , aa = 0; end
    bb = theta_new+delta; if bb > 2*pi , bb = 2*pi; end
    ttt = linspace(aa, bb, 101);
    P = create_ellipse_proj(R, ttt, T_pe_H);
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
    
    theta_new = ttt(min_idx);
    tt_new = t(min_idx);
    local_min_point = [ppx(min_idx), ppy(min_idx), ppz(min_idx)].'
    
    k=k+1;
end
lowest_point = local_min_point;
plot3(lowest_point(1,:), lowest_point(2,:), lowest_point(3,:), 'b*');