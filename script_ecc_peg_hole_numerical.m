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
%   2015-09-25: Fixed the bug when both theta_x and theta_y are small, the
%               equation will fail to be solved.

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
theta_x_pe_H = 0.0;       % ���ƿ�(H)����ϵx��ת���ĽǶ�
theta_y_pe_H = 0.10472;         % ���ƿ�(H)����ϵy��ת���ĽǶ�
theta_z_pe_H = 0;                % ���ƿ�(H)����ϵz��ת���ĽǶ�(���趨����Ϊ����Ҫ�ڸ��᷽����ת)
O_x_pe_H = -1.5;                         % ��ĵ������ĵ��ڿ�����ϵx�᷽��ƫ��׵����ĵ�ľ���
O_y_pe_H = 0.0;                         % ��ĵ������ĵ��ڿ�����ϵy�᷽��ƫ��׵����ĵ�ľ���
O_z_pe_H = 5;                         % ��ĵ������ĵ��ڿ�����ϵz�᷽��ƫ��׵����ĵ�ľ���

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
theta = linspace(0, 2*pi, 721); % ���̫��Ӱ�����Ӵ�����жϾ���
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
% [  cos(ty), sin(tx)*sin(ty), cos(tx)*sin(ty), px]
% [        0,         cos(tx),        -sin(tx), py]
% [ -sin(ty), cos(ty)*sin(tx), cos(tx)*cos(ty), pz]
% [        0,               0,               0,  1]
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

%% ����͵�

% legacy code-----------------------------------------------------------------
% % ����ĵ���Բ��Z_h�����ϵ���͵�
% [~, j]=find(points == min(min(points(3,:))));

% 20140716-----------------------------------------------------------------
axis_a_pe_H_proj = [axis_a_pe_H(1:2); 0]; % z����XoYƽ���ͶӰ
costt = (axis_a_pe_H_proj' * axis_a_pe_H)/(norm(axis_a_pe_H_proj)*norm(axis_a_pe_H)); %z����XoYƽ��ļн�

vector_R = [0 0 0]'; % �����ĵ㵽��͵������
vector_R(3) = -R_hole * costt; % z���������������
a_R = axis_a_pe_H_proj / norm(axis_a_pe_H_proj) * R_hole * sqrt(1-costt^2);
vector_R(1:2) = a_R(1:2); % x y������ͶӰ�������

lowest_point = O_pe_H + vector_R;

% ����ĵ���Բ��Z_h�����ϵ���͵�
plot3([lowest_point(1) lowest_point(1)],...
    [lowest_point(2) lowest_point(2)],...
    [lowest_point(3) 0], 'r');
hold on
plot3([lowest_point(1) lowest_point(1)],...
    [lowest_point(2) lowest_point(2)],...
    [lowest_point(3) 0], 'ro');
%% �ų���͵�ͶӰ����hole�ڵ����
if norm(lowest_point(1:2)) > R_hole
    disp('lowest point OUT of range!');
    d_lowest = NaN;
    flag.outrange = 1;
    return;
end

%%
% �����Բ��ͶӰԲ�Ľ���
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
% % f1 = x^2 + y^2 - R^2;
% A1 = 1; B1 = 1; C1 = 0; D1 = 0; E1 = 0; F1 = -R_hole^2;
% % f2 = (x - O_x_pe_H)^2*Py^2 - 2*(x - O_x_pe_H)*(y - O_y_pe_H)*Q + (y - O_y_pe_H)^2*(Q^2+Px^2) - R^2;
% %A2 = Py^2;
% %B2 = Px^2 + Q^2;
% %C2 = - Q;
% %D2 = Q*O_y_pe_H - Py^2*O_x_pe_H;
% %E2 = Q*O_x_pe_H - (Q^2 + Px^2)*O_y_pe_H;
% %F2 = Py^2*O_x_pe_H^2 + Px^2*O_y_pe_H^2 + Q^2*O_y_pe_H^2 - 2*Q*O_x_pe_H*O_y_pe_H - R^2;
%
% % ֱ��ʹ�÷�������������Ϊ����Բ��ȡ�㣬Ȼ����ϵ���ķ���
% sel_points = [1 73 145 217 289];
% coe_xx = points(1,sel_points)';
% coe_yy = points(2,sel_points)';
%
% coe_C = [coe_xx.^2, coe_yy.^2, 2*coe_xx.*coe_yy, 2*coe_xx, 2*coe_yy];
%
% % ��ϵ���������ȣ������ͬ���칹
% if rank(coe_C) < 5
%     coe_C = [coe_C ; zeros(1, size(coe_C, 2))];
%     sol_C = coe_C \ ones(size(coe_xx)+1);
% else
%     sol_C = coe_C\ones(size(coe_xx)); % ͨ����Բ�ϵ��������A��B��C��D��E���ϵ����F=1
% end
%
% A2 = sol_C(1);
% B2 = sol_C(2);
% C2 = sol_C(3);
% D2 = sol_C(4);
% E2 = sol_C(5);
% F2 = -1;
% f1 = [A1 C1 D1; C1 B1 E1; D1 E1 F1];
% f2 = [A2 C2 D2; C2 B2 E2; D2 E2 F2];
% S = intersectConics(f1, f2);
% % plot(S(1,:) ./ S(3,:) , S(2,:) ./ S(3,:), 'ro');
% [m_s, ~] = size(S);
% if m_s == 3
%     x = S(1,:) ./ S(3,:);
%     y = S(2,:) ./ S(3,:);
% else
%     x=[];
%     y=[];
% end

% 2015-09-28 ��ֵ������Բ����Բ����
P = create_ellipse_proj(R, theta, T_pe_H); % ����ת���peg����㼯
xe = P(1,:);
ye = P(2,:);
f = abs(xe.^2 + ye.^2 - R^2); % ����Բ�ϵĵ����Բ�ķ��̣���Сֵ��

local_min_idx = localMaximum(-f);
local_min_n = length(local_min_idx);
local_min = theta(local_min_idx);

delta = pi/18; % 10 degree
n_step = 100;
k=1;
while(k < 10) % ����10�λ�����������ƽ�
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
        g_idx(i) = t_idx(1); % ��ֹt_idx��ͬһλ�����������
        
        local_min(i) = tt(g_idx(i));
    end
    k = k+1;
end

P = create_ellipse_proj(R, local_min, T_pe_H);
xe = P(1,:);
ye = P(2,:);

f = abs(xe.^2 + ye.^2 - R^2);

% ɾȥ���ǽ�����е�ĵ�
f_idx = f > 10^-6;
xe(f_idx)=[];
ye(f_idx)=[];
f(f_idx) = [];

x = xe;
y = ye;

%% �жϽ������
% ͶӰԲ��hole�Ľ�����0����1����2����3����4��5���������Ӧ��
% �Ӵ��������¼��������
% 	1. �޽Ӵ���                        (0������)
% 	2. ��һ���Ӵ���
% 	2.1  ��ĵ���Բ�������             (1������/3������)
% 	2.2  ��ĵ���Բ��ײ�����           (2������)
%   2.3  ��Ĳ��������һ���Ӵ���       (2������)
% 	3. �������Ӵ���                    (2������/4������)
% 	4. �������Ӵ���                    (2������/4������)

solution_count = length(x);
flag.inter = solution_count;
if solution_count > 4 % ����Բ�غϣ��������insertion
    d_lowest = -h_ecc;
    flag.inter = Inf;
    return;
end

% ���ÿ��ͶӰ�����Ӧ��peg�����ϵ�ԭ��ĸ߶�z��ͨ���Ƚϸ߶��ж�ʵ�ʽӴ����м���
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
        if norm(z_cmp) > 10^-6 % ��������һ��һ�� - ֻ��һ���Ӵ���
            z = z(idx(1));
            x = x(idx(1));
            y = y(idx(1));
            flag.contact = 1;
        else
            flag.contact = 2;
        end
    case 3
        if norm(z_cmp(1:2)) < 10^-6 % ��������һ���� - �������Ӵ���
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
        if norm(z_cmp) < 10^-6 % �ĸ�����һ���� - ���ڿ���װ����ȥ
            d_lowest = 0;
            flag.contact = 4;
            return;
        elseif norm(z_cmp(1:2)) < 10^-6 % ��������һ���� - �������Ӵ���
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

%% ���������Ƿ��нӴ���

% peg����Բ�ϵĵ�
px = points(1,:);
py = points(2,:);
pz = points(3,:);

% ��peg�����������᷽��ļн�cos(alpha) cos(beta) cos(gamma)
ca = axis_o_pe' * axis_a_pe_H;
cb = axis_n_pe' * axis_a_pe_H;
cc = axis_a_pe' * axis_a_pe_H;
cs = ca^2 + cb^2;

% peg�����������
% x = x' + t*cos(alpha)
% y = y' + t*cos(beta)
% z = z' + t*cos(gamma)
% t�ǲ�����(x',y',z')������������һ��Բ�����ϵĵ㣬����(x',y',z')��peg������ȡ��
cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);

% ȥ�����б�ʽС����ĵ�Ͷ�Ӧ��ֵ
delta_idx = cdelta<0;
px(delta_idx) = NaN;
py(delta_idx) = NaN;
pz(delta_idx) = NaN;
cdelta(delta_idx) = NaN;

% �����Ӧ�Ĳ���t��ȡֵ.��ʱ���ϲ�������ʵ��������Բ���Ľ��ߵĲ�������
t = [(-(px*ca+py*cb) + sqrt(cdelta)) / cs...
    (-(px*ca+py*cb) - sqrt(cdelta)) / cs];

% ��������ϵĵ�.������t<0ʱ��Ӧ����peg�������µĽ��ߣ�ʵ������в����ڣ��ʰ���
% Ӧ�ĵ�Ͷ�Ӧ��ֵ��ȥ
ppx = [px px] + t * ca;
ppy = [py py] + t * cb;
ppz = [pz pz] + t * cc;

t_idx = t<0;
ppx(t_idx) = NaN;
ppy(t_idx) = NaN;
ppz(t_idx) = NaN;

[~, idx] = min(ppz);
if idx <= length(t)/2
    delta_flag = 1; % ����t�Ľ�ķ��ŵı�־λ
else
    delta_flag = 2;
end

thetatheta = [theta theta]; % t�Ľ�����֧,��Ӧԭ��������Ҫ����
theta_new = thetatheta(idx);
theta_range = zeros(1, 101);

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

% �Ƚ�һ�²����͵����ĸ߶�
edge_x = edge_min_point(1);
edge_y = edge_min_point(2);
edge_z = edge_min_point(3);

contact_point = [x;y;z];

if abs(edge_z - z(1)) > 10^-6 % �����͵����߶Ȳ�ͬ
    if edge_z < z(1) % ����Ӵ�
        h_down = edge_z;
        contact_point = edge_min_point;
        flag.contact = -1;
    else % ����Ӵ�
        h_down = z(1);
    end
else % ����͵���߶���ͬ
    if flag.contact == 2 && abs(t(min_idx)) > 10^-6 % ����Ӵ�,�������Ͳ�����ͬһ����
        contact_point = [contact_point edge_min_point];
        flag.contact = 3;               
    end
    h_down = z(1);
end

d_lowest = lowest_point(3) - h_down;
plot3(contact_point(1,:), contact_point(2,:), contact_point(3,:), 'r*');

%% �쳣����
% if d_lowest > 5 || d_lowest < -5
%     d_lowest = NaN;
%     fprintf(2,'Unexpected result!\n');
% end

%% ��͵���ֵ
title(['d_{lowest} = ' num2str(d_lowest)]);
disp(flag);





