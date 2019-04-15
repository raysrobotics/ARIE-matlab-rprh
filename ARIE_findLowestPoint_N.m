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
%   Rui Li (raysworld@outlook.com)
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

% ���λ
flag = [];
flag.outrange = 0;          % ��͵��Ƿ���hole��Χ�ڣ�����Ϊ0��������Ϊ1
flag.inter = 0;             % ͶӰ���Բ���������ֵΪn����n������
flag.contact = 0;           % ����׵ĽӴ��������ֵΪn����n���Ӵ��㣬-1Ϊ����Ӵ�
% �����͵��z����ֵ
% d_lowest = NaN;
% ƽ����
% x��yΪ�仯����zΪ������ڿ�����ƽ��̧��߶ȣ���Ϊһ����ֵ���ɣ�����5
O_x_pe_H = param(1);        % x
O_y_pe_H = param(2);        % y
O_z_pe_H = param(3);        % z

%% ����
% ��װ뾶
R = options.R;
% ƫ�ľ�(��X_pe�᷽��)
% d_ecc = options.d_ecc;
% ƫ�ľ�(��Z_pe����)
h_ecc = options.h_ecc;

% peg-e��������hole���������ɼн�
theta_x_pe_H = options.theta_x;
theta_y_pe_H = options.theta_y;
% theta_z_pe_H = options.theta_z;

%% ��׼Բ�����ڸ��ֱ任�ͻ�ͼ
theta = linspace(0, 2*pi, 721);
rr = R;
xx = rr * cos(theta);
yy = rr * sin(theta);
zz = zeros(1,length(theta));
ww = ones(1,length(theta)); % �������

%%
% hole����ϵ(������ϵ)------------------------------------------------------
% frame_H = eye(4);
if display == 2 || display == 3
    % ����hole
    plot3(xx, yy, zz, 'k'); hold on
    % �������ע
    xlabel('x'); ylabel('y'); zlabel('z');

    grid on; axis equal
    axis([-2*R 2*R -2*R 2*R -2 8]);
end

% peg-e����ϵ(ƫ�����²��̶�)-----------------------------------------------
% peg-e��hole����ϵ�Ĺ�ϵ����[Xpe Ype Zpe theta_xe theta_ye theta_ze]
axis_o_pe = [1 0 0]';
axis_n_pe = [0 1 0]';
axis_a_pe = [0 0 1]';
O_pe = [0 0 0]';
frame_pe = [axis_o_pe axis_n_pe axis_a_pe O_pe;0 0 0 1];
% peg-e����ϵ��hole����ϵ�µı�ʾ(�任)
T_pe_H = trans([O_x_pe_H O_y_pe_H O_z_pe_H])*rot([0 1 0], theta_y_pe_H)*rot([1 0 0], theta_x_pe_H);
frame_pe_H = T_pe_H * frame_pe;
% peg-e����ϵ�������ᡢԭ����hole����ϵ�µı�ʾ(�任���)
% axis_o_pe_H = frame_pe_H(1:3,1);
% axis_n_pe_H = frame_pe_H(1:3,2);
axis_a_pe_H = frame_pe_H(1:3,3);
O_pe_H = frame_pe_H(1:3,4);

%% ��peg-e�����Բ
% peg-e�����Բ
points = T_pe_H * [xx;yy;zz;ww];
if display == 2 || display == 3
    plot3(points(1,:),points(2,:),points(3,:),'b');
    hold on
    % peg-e�����Բ��ԭ��
    plot3(O_pe_H(1),O_pe_H(2),O_pe_H(3),'r.');
    % peg-e�����Բ��XOY(hole)ƽ���ϵ�ͶӰ
    plot3(points(1,:),points(2,:),zz,'b');
end
%% ����͵�

% 20140716-----------------------------------------------------------------
axis_a_pe_H_proj = [axis_a_pe_H(1:2); 0]; % z����XoYƽ���ͶӰ
costt = (axis_a_pe_H_proj' * axis_a_pe_H)/(norm(axis_a_pe_H_proj)*norm(axis_a_pe_H)); %z����XoYƽ��ļн�

vector_R = [0 0 0]'; % �����ĵ㵽��͵������
vector_R(3) =  -R * costt; % z���������������
a_R = axis_a_pe_H_proj / norm(axis_a_pe_H_proj) * R * sqrt(1-costt^2);
vector_R(1:2) = a_R(1:2); % x y������ͶӰ�������

lowest_point = O_pe_H + vector_R;

if display == 2 || display == 3
    % ����ĵ���Բ��Z_h�����ϵ���͵�
    plot3([lowest_point(1) lowest_point(1)],...
          [lowest_point(2) lowest_point(2)],...
          [lowest_point(3) 0], 'r');
    hold on
    plot3([lowest_point(1) lowest_point(1)],...
          [lowest_point(2) lowest_point(2)],...
          [lowest_point(3) 0], 'ro');
end

%% �ų���͵�ͶӰ����hole�ڵ����
if norm(lowest_point(1:2)) > R
	if (display == 1 || display == 3), disp('lowest point OUT of range!'), end
    d_lowest = NaN;
    O_pe_H(3) = O_pe_H(3) - lowest_point(3);
    
    flag.outrange = 1;
    return;
end

%% ��Բ����Բ����

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

% ��β����Ŀ�����ͬһ����
if (xe(1)-xe(end))^2+(ye(1)-ye(end))^2 < 10^-6
    xe(end)=[];
    ye(end)=[];
end

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
% ppx = [px px] + t * ca;
% ppy = [py py] + t * cb;
ppz = [pz pz] + t * cc;

t_idx = t<0;
% ppx(t_idx) = NaN;
% ppy(t_idx) = NaN;
ppz(t_idx) = NaN;

[~, idx] = min(ppz);
if idx <= length(t)/2
    delta_flag = 1; % ����t�Ľ�ķ��ŵı�־λ
else
    delta_flag = 2;
end

thetatheta = [theta theta]; % t�Ľ�����֧,��Ӧԭ��������Ҫ����
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

% �Ƚ�һ�²����͵����ĸ߶�
% edge_x = edge_min_point(1);
% edge_y = edge_min_point(2);
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
O_pe_H(3) = O_pe_H(3) - h_down;

%% Add title to the figure
if display == 2 || display == 3
    plot3(contact_point(1,:), contact_point(2,:), contact_point(3,:), 'r*');
    title(['d_{lowest} = ' num2str(d_lowest)]);
    disp(flag);
end

