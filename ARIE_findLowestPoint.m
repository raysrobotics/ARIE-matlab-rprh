function [flag, d_lowest, O_pe_H] = ARIE_findLowestPoint(param, options, display)
%ARIE_findLowestPoint   Find lowest point of an Attractive Region in 
%                       Environment.
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
%   2015-09-25: Fixed the bug when both theta_x and theta_y are small, the
%               equation will fail to be solved.

% ���λ
flag = [];
flag.outrange = 0;          % ��͵��Ƿ���hole��Χ�ڣ�����Ϊ0��������Ϊ1
flag.inter = 0;             % ͶӰ���Բ���������ֵΪn����n������
flag.contact = 0;           % ����׵ĽӴ��������ֵΪn����n���Ӵ��㣬-1Ϊ����Ӵ�
% �����͵��z����ֵ
d_lowest = NaN;
% ƽ����
% x��yΪ�仯����zΪ������ڿ�����ƽ��̧��߶ȣ���Ϊһ����ֵ���ɣ�����5
O_x_pe_H = param(1);        % x
O_y_pe_H = param(2);        % y
O_z_pe_H = param(3);        % z

%% ����
% ��װ뾶
R = options.R;
% ƫ�ľ�(��X_pe�᷽��)
d_ecc = options.d_ecc;
% ƫ�ľ�(��Z_pe����)
h_ecc = options.h_ecc;

% peg-e��������hole���������ɼн�
theta_x_pe_H = options.theta_x;
theta_y_pe_H = options.theta_y;
theta_z_pe_H = options.theta_z;

%% ��׼Բ�����ڸ��ֱ任�ͻ�ͼ
t = linspace(0, 2*pi, 360);
rr = R;
xx = rr * cos(t);
yy = rr * sin(t);
zz = zeros(1,length(t));
ww = ones(1,length(t)); % �������

%%
% hole����ϵ(������ϵ)------------------------------------------------------
frame_H = eye(4);
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
axis_o_pe_H = frame_pe_H(1:3,1);
axis_n_pe_H = frame_pe_H(1:3,2);
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
% legacy code--------------------------------------------------------------
% % ����ĵ���Բ��Z_h�����ϵ���͵�
% [~, j]=find(points == min(min(points(3,:))));
% if display == 2 || display == 3
%     % ����ĵ���Բ��Z_h�����ϵ���͵�
%     plot3([points(1,j) points(1,j)],...
%         [points(2,j) points(2,j)],...
%         [points(3,j) zeros(1,length(j))],...
%         'r');
% end
% lowest_point = [points(1,j(1)) points(2,j(1)) points(3,j(1))].';

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
    O_pe_H(3) = NaN;
    
    flag.outrange = 1;
    return;
end

%% �󽻵�
% legacy code--------------------------------------------------------------
% syms x y
% Px = sin(theta_x_pe_H)*tan(theta_x_pe_H) + cos(theta_x_pe_H);
% Py = sin(theta_y_pe_H)*tan(theta_y_pe_H) + cos(theta_y_pe_H);
% Q = tan(theta_x_pe_H)*tan(theta_y_pe_H);

% f1 = x^2 + y^2 - R^2;
% f2 = (x - O_x_pe_H)^2*Py^2 - 2*(x - O_x_pe_H)*(y - O_y_pe_H)*Q + (y - O_y_pe_H)^2*(Q^2+Px^2) - R^2;
% 
% S = solve(f1,f2, 'Real', true);
% 
% A = (S.x == real(S.x));
% x = eval(S.x(A));
% B = (S.y == real(S.y));
% y = eval(S.y(B));

% 20140509-----------------------------------------------------------------
% % f1 = x^2 + y^2 - R^2;
% A1 = 1; B1 = 1; C1 = 0; D1 = 0; E1 = 0; F1 = -R^2;
% % f2 = (x - O_x_pe_H)^2*Py^2 - 2*(x - O_x_pe_H)*(y - O_y_pe_H)*Q + (y - O_y_pe_H)^2*(Q^2+Px^2) - R^2;
% % A2 = Py^2;
% % B2 = Px^2 + Q^2;
% % C2 = - Q;
% % D2 = Q*O_y_pe_H - Py^2*O_x_pe_H;
% % E2 = Q*O_x_pe_H - (Q^2 + Px^2)*O_y_pe_H;
% % F2 = Py^2*O_x_pe_H^2 + Px^2*O_y_pe_H^2 + Q^2*O_y_pe_H^2 - 2*Q*O_x_pe_H*O_y_pe_H - R^2;
% 
% % ֱ��ʹ�÷�������������Ϊ����Բ��ȡ�㣬Ȼ����ϵ���ķ���
% sel_points = [1 61 121 181 241];
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
% 
% f1 = [A1 C1 D1; C1 B1 E1; D1 E1 F1];
% f2 = [A2 C2 D2; C2 B2 E2; D2 E2 F2];
% S = intersectConics(f1, f2);
% % plot(S(1,:) ./ S(3,:) , S(2,:) ./ S(3,:), 'ro');
% [m_s, ~] = size(S);
% if m_s == 3
%     x = S(1,:) ./ (S(3,:)+eps);
%     y = S(2,:) ./ (S(3,:)+eps);
% else
%     x=[];
%     y=[];
% end

% 2015-09-28 ��ֵ������Բ����Բ����
P = create_ellipse_proj(R, t, T_pe_H); % ����ת���peg����㼯
xe = P(1,:);
ye = P(2,:);
f = abs(xe.^2 + ye.^2 - R^2); % ����Բ�ϵĵ����Բ�ķ��̣���Сֵ��

local_min_idx = localMaximum(-f);
local_min_n = length(local_min_idx);
local_min = t(local_min_idx);

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
% �Ƚ��ĸ�������������Ϊ��������
if solution_count == 4 % �Ƚ��ĸ���ĸ߶ȣ����¸߶ȵ͵�����
    if (display == 1 || display == 3), disp('FOUR intersection points!'), end
    flag.inter = 4;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    [~, idx] = sort(z);
    z = z(idx(1:2));
    x = x(idx(1:2));
    y = y(idx(1:2));
end
solution_count = length(x); % 2015-09-28

% if solution_count == 2 % �Ƚ��������Ƿ���ͬһ����(2014-07-18)
%     z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
%     if norm([x(1)-x(2) y(1)-y(2) z(1)-z(2)]) < 10^-6
%         x = x(1);
%         y = y(1);
%         z = z(1);
%     end
% end
% solution_count = length(x); % 2015-09-28

% 0������
if solution_count == 0
    if (display == 1 || display == 3), disp('NO intersection points!'), end
    flag.inter = 0;
    
elseif solution_count == 1 
    if (display == 1 || display == 3), disp('ONE intersection point!'), end
    flag.inter = 1;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    h_down = z;
    d_lowest = lowest_point(3) - h_down;
% 20140718-----------------------------------------------------------------
    if (display == 1 || display == 3), disp('ONE contact point!'), end
    flag.contact = 1;
% 20140718-END-------------------------------------------------------------
elseif solution_count == 3
    if (display == 1 || display == 3), disp('THREE intersection points!'), end
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
        if (display == 1 || display == 3), disp('TWO intersection points!'), end
        flag.inter = 2;
    end
    
    % ���㽻��λ�õ�holeƽ��ľ���
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    if display == 2 || display == 3
        plot3([x(1) x(1)],[y(1) y(1)],[0 z(1)], 'g-');
        plot3([x(2) x(2)],[y(2) y(2)],[0 z(2)], 'g-');
    end
    
    % ���潻�����굽�µı�����
    inter_x = x;
    inter_y = y;
    inter_z = z;
    
    % ���Ӵ�������Ӵ���
    if norm(z(1)-z(2)) < 10^-4 % ����߶�����˵�����ܴ��ڵ������Ӵ��㣬��Ҫ��һ���ж�
        h_down = z(1);
        % �����ж��Ƿ���ڵ������Ӵ���
        % ����һ����ֱ��XoYƽ���ƽ��,��ƽ���Ope�����͵�
%         % ƽ�淽��ΪAx+By=1����
%         % coe_planar(1)*x + coe_planar(2)*y - coe_planar(3)*1 = 0
%         % �������� - O_pe_H��lowest_point���룬���ɽ��A��B
%         if abs(O_pe_H(1)-lowest_point(1))<10^-6%0.02  % ƽ��ӽ���X���y���غ�ʱ�����̲�̬���������棬��ֱ�Ӽ���
%             coe_planar = [1,0,O_pe_H(1)];
%         elseif abs(O_pe_H(2)-lowest_point(2))<10^-6%0.02
%             coe_planar = [0,1,O_pe_H(2)];
%         else
%             coe_planar = [O_pe_H(1) O_pe_H(2);lowest_point(1) lowest_point(2)]\[1;1];
%             coe_planar(3) = 1;
%         end
        % 2015-09-28
        % ƽ�淽��Ϊ Ax+By+C=0
        % A = Y2 - Y1
        % B = X1 - X2
        % C = X2*Y1 - X1*Y2
        p1 = O_pe_H(1:2);
        p2 = lowest_point(1:2);
        A = p2(2) - p1(2);
        B = p1(1) - p2(1);
        C = p2(1)*p1(2) - p1(1)*p2(2);
        
        coe_planar = [A B C];
        
        cy = cos(theta_y_pe_H);
        sy = sin(theta_y_pe_H);
        cx = cos(theta_x_pe_H);
        sx = sin(theta_x_pe_H);
        syms x y
            f11  = (x - O_x_pe_H)*cy-(h_down-O_z_pe_H)*sy;
            f121 = (x - O_x_pe_H)*sx*sy;
            f122 = (y - O_y_pe_H)*cx;
            f123 = (h_down-O_z_pe_H)*sx*cy;
        f1 = f11^2 + (f121 + f122 + f123)^2 - R^2;        
        f2 = coe_planar(1)*x + coe_planar(2)*y + coe_planar(3);
        
        S = solve(f2, f1, 'Real', 1);       
        
        % 2015-09-25 S.x & S.y   
        if isempty(S)
            disp('Equation solve failed!');
            d_lowest = NaN;
            flag.outrange = 1;
            return;
        elseif isempty(S.x) || isempty(S.y)
            d_lowest = -O_z_pe_H;
            flag.inserted = 1;
            return; 
        end
        
        A = (S.x == real(S.x));
        x = eval(S.x(A));
        B = (S.y == real(S.y));
        y = eval(S.y(B));
        
%         plot3(x(1), y(1), h_down, 'go');
%         plot3(x(2), y(2), h_down, 'go');
        
        % �ж��ĸ��������ϵĵ�(Ope����͵㡢���ϵ��������߳�ֱ��)
        v1 = O_pe_H - lowest_point;
        v2 = [x(1);y(1);h_down] - lowest_point;
        if abs(v1.' * v2) < 10^-6
            edge_point = [x(1) y(1) h_down].';
        else
            edge_point = [x(2) y(2) h_down].';
        end
        if (display == 2 || display == 3)
            plot3(edge_point(1), edge_point(2), edge_point(3), 'go');
        end
        
        % �ж�����һ���ͶӰ��hole��λ�ù�ϵ
        if abs(edge_point(1)^2 + edge_point(2)^2 - R^2) <= 10^-4 % ����hole�ϣ�����Ӵ�
            if (display == 1 || display == 3), disp('three-point contact!'), end
            flag.contact = 3;
            d_lowest = lowest_point(3) - h_down;
        elseif edge_point(1)^2 + edge_point(2)^2 - R^2 < 0 % ����hole�ڣ�����Ӵ�
            if (display == 1 || display == 3), disp('two-point contact!'), end
            flag.contact = 2;
            d_lowest = lowest_point(3) - h_down;
        elseif edge_point(1)^2 + edge_point(2)^2 - R^2 > 0 % ����hole�⣬�����Ӵ�
            if (display == 1 || display == 3), disp('edge-only contact!'), end
            flag.contact = -1;
            d1 = norm(lowest_point(1:2) - edge_point(1:2));
            mid_point = [sum(inter_x)/2  sum(inter_y)/2  sum(inter_z)/2].';
            
            d2_cos = ([0;0] - mid_point(1:2))' * (edge_point(1:2) - mid_point(1:2));
            if d2_cos >= 0 % mid_point��O_hole��edge_point_proj��ͬ��
                d2 = norm(mid_point(1:2) - edge_point(1:2)) - norm(mid_point(1:2)) - R;
            else % mid_point��O_hole��edge_point_proj���м�
                d2 = norm(mid_point(1:2) - edge_point(1:2)) + norm(mid_point(1:2)) - R;
            end
            
            d3 = abs(h_down - lowest_point(3));
            d4 = d2*d3/d1;
            
            d_lowest = lowest_point(3) - (h_down - d4);
        end
    else % ����߶����ϴ�˵��ֻ��һ���Ӵ��㣬�˴����ж�һ���ǲ��ǲ���Ӵ�(2014-08-22)        
        % ȡ��peg����Բ�ϵĵ�
        px = points(1,:);
        py = points(2,:);
        pz = points(3,:);
        
        % peg��z����hole����ϵ���ɵķ�������
        ca = axis_o_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cb = axis_n_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cc = axis_a_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cs = ca^2 + cb^2;
        
        % �Ƶ�����delta�����������Բ�����潻���ϵĵ�
        cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);

        idx = cdelta>=0;
        
        px = px(idx);
        py = py(idx);
        pz = pz(idx);
        cdelta = cdelta(idx);
        
        % t1Ϊ�������ߵ�����һ֧
        t1 = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;
        
        % [ppx, ppy, ppz]����Բ���潻���ϵĵ�
        ppx = px + t1 * ca;
        ppy = py + t1 * cb;
        ppz = pz + t1 * cc;

        if abs(O_pe_H(1)-lowest_point(1))<10^-6%0.02  % ƽ��ӽ���X���y���غ�ʱ�����̲�̬���������棬��ֱ�Ӽ���
            coe_planar = [1,0,O_pe_H(1)];
        elseif abs(O_pe_H(2)-lowest_point(2))<10^-6%0.02
            coe_planar = [0,1,O_pe_H(2)];
        else
            coe_planar = [O_pe_H(1) O_pe_H(2);lowest_point(1) lowest_point(2)]\[1;1];
            coe_planar(3) = 1;
        end
        
        cd = abs(coe_planar(1)*ppx + coe_planar(2)*ppy - coe_planar(3)*1) < 10^-4;
        [~, idx] = min(ppz(cd));  % [ppx(idx), ppy(idx), ppz(idx)]�����ܵĲ���Ӵ���

        if isempty(idx) || (min(z) < ppz(idx)) % һ��Ӵ�            
            if (display == 1 || display == 3), disp('one-point contact!'), end
            flag.contact = 1;

            h_down = min(z);
        else % �����Ӵ�
            if (display == 1 || display == 3), disp('edge-only contact!'), end
            flag.contact = -1;
            
            h_down = ppz(idx);            
        end
        d_lowest = lowest_point(3) - h_down;
        
        if d_lowest>0
            d_lowest = 0;
        end
    end
end

% calculate the z-position of O_pe_H
O_pe_H(3) = O_pe_H(3) - h_down;

%% Add title to the figure
if display == 2 || display == 3
    title(['d_{lowest} = ' num2str(d_lowest)]);
    disp(flag);
end

