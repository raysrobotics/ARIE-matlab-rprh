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
%   2014-07-18: Distinguish radius of peg and radius of hole
%   2014-08-22: Improvements on edge-only detect
%

%% 
clear all
close all
clc
%% ����
% ��ײ���
R = 30;      % ��װ뾶
R_peg = 30.0;
R_hole = 30.2;
d_ecc = 1;  % ƫ�ľ�(��X_pe�᷽��)
h_ecc = 3;  % ƫ�ľ�(��Z_pe����)

% peg-e��hole����ϵ�Ĺ�ϵ����[Xpe Ype Zpe theta_xe theta_ye theta_ze]
theta_x_pe_H = pi/3;       % ���ƿ�(H)����ϵx��ת���ĽǶ�
theta_y_pe_H = 0;         % ���ƿ�(H)����ϵy��ת���ĽǶ�
theta_z_pe_H = 0;                % ���ƿ�(H)����ϵz��ת���ĽǶ�(���趨����Ϊ����Ҫ�ڸ��᷽����ת)
O_x_pe_H = -0.4667;                         % ��ĵ������ĵ��ڿ�����ϵx�᷽��ƫ��׵����ĵ�ľ���
O_y_pe_H = 0.125;                         % ��ĵ������ĵ��ڿ�����ϵy�᷽��ƫ��׵����ĵ�ľ���
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
t = linspace(0, 2*pi, 360);
% rr = R;
% xx = rr * cos(t);
% yy = rr * sin(t);
% zz = zeros(1,length(t));
% ww = ones(1,length(t));

rr_p = R_peg;
xx_p = rr_p * cos(t);
yy_p = rr_p * sin(t);
zz_p = zeros(1,length(t));
ww_p = ones(1,length(t));

rr_h = R_hole;
xx_h = rr_p * cos(t);
yy_h = rr_p * sin(t);
zz_h = zeros(1,length(t));
ww_h = ones(1,length(t));


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

%% ����͵�

% legacy code-----------------------------------------------------------------
% % ����ĵ���Բ��Z_h�����ϵ���͵�
% [~, j]=find(points == min(min(points(3,:))));

% 20140716-----------------------------------------------------------------
axis_a_pe_H_proj = [axis_a_pe_H(1:2); 0]; % z����XoYƽ���ͶӰ
costt = (axis_a_pe_H_proj' * axis_a_pe_H)/(norm(axis_a_pe_H_proj)*norm(axis_a_pe_H)); %z����XoYƽ��ļн�

vector_R = [0 0 0]'; % �����ĵ㵽��͵������
vector_R(3) =  -R_hole * costt; % z���������������
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
% f1 = x^2 + y^2 - R^2;
A1 = 1; B1 = 1; C1 = 0; D1 = 0; E1 = 0; F1 = -R_hole^2;
% f2 = (x - O_x_pe_H)^2*Py^2 - 2*(x - O_x_pe_H)*(y - O_y_pe_H)*Q + (y - O_y_pe_H)^2*(Q^2+Px^2) - R^2;
%A2 = Py^2;
%B2 = Px^2 + Q^2;
%C2 = - Q;
%D2 = Q*O_y_pe_H - Py^2*O_x_pe_H;
%E2 = Q*O_x_pe_H - (Q^2 + Px^2)*O_y_pe_H;
%F2 = Py^2*O_x_pe_H^2 + Px^2*O_y_pe_H^2 + Q^2*O_y_pe_H^2 - 2*Q*O_x_pe_H*O_y_pe_H - R^2;

% ֱ��ʹ�÷�������������Ϊ����Բ��ȡ�㣬Ȼ����ϵ���ķ���
sel_points = [1 73 145 217 289];
coe_xx = points(1,sel_points)';
coe_yy = points(2,sel_points)';

coe_C = [coe_xx.^2, coe_yy.^2, 2*coe_xx.*coe_yy, 2*coe_xx, 2*coe_yy];
sol_C = coe_C\ones(size(coe_xx)); % ͨ����Բ�ϵ��������A��B��C��D��E���ϵ����F=1

A2 = sol_C(1);
B2 = sol_C(2);
C2 = sol_C(3);
D2 = sol_C(4);
E2 = sol_C(5);
F2 = -1;
f1 = [A1 C1 D1; C1 B1 E1; D1 E1 F1];
f2 = [A2 C2 D2; C2 B2 E2; D2 E2 F2];
S = intersectConics(f1, f2);
% plot(S(1,:) ./ S(3,:) , S(2,:) ./ S(3,:), 'ro');
[m_s, ~] = size(S);
if m_s == 3
    x = S(1,:) ./ S(3,:);
    y = S(2,:) ./ S(3,:);
else
    x=[];
    y=[];
end

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
    disp('four intersection points!');
	flag.inter = 4;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    [~, idx] = sort(z);
    z = z(idx(1:2));
    x = x(idx(1:2));
    y = y(idx(1:2));
end
if solution_count == 2 % �Ƚ��������Ƿ���ͬһ����(2014-07-18)
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    if norm([x(1)-x(2) y(1)-y(2) z(1)-z(2)]) < 10^-6
        x = x(1);
        y = y(1);
        z = z(1);
    end
end

solution_count = length(x);
% 0������
if solution_count == 0
    disp('NO intersection points!');
    flag.inter = 0;
elseif solution_count == 1
    disp('ONE intersection point!');
	flag.inter = 1;
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    
    h_down = z;
    d_lowest = lowest_point(3) - h_down;
    
% 20140718-----------------------------------------------------------------
    disp('ONE contact point!');
    flag.contact = 1;
% 20140718-END-------------------------------------------------------------
elseif solution_count == 3
    disp('THREE intersection points!');
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
        disp('TWO intersection points!');
        flag.inter = 2;
    end
    
    % ���㽻��λ�õ�holeƽ��ľ���
    z = ((y - O_y_pe_H)*sin(theta_x_pe_H) - (x - O_x_pe_H)*cos(theta_x_pe_H)*sin(theta_y_pe_H))/(cos(theta_x_pe_H)*cos(theta_y_pe_H))+O_z_pe_H;
    plot3([x(1) x(1)],[y(1) y(1)],[0 z(1)], 'g-');
    plot3([x(2) x(2)],[y(2) y(2)],[0 z(2)], 'g-');
    
    % ���潻�����굽�µı�����
    inter_x = x;
    inter_y = y;
    inter_z = z;
    
    if norm(z(1)-z(2)) < 10^-4 % ���Ӵ�������Ӵ���
        h_down = z(1);
        % �����ж��Ƿ���ڵ������Ӵ���
        % ����һ����ֱ��XoYƽ���ƽ��,��ƽ���Ope�����͵�
        % ƽ�淽��ΪAx+By=1
        if abs(O_pe_H(1)-lowest_point(1))<0.02
            coe_planar = [1,0,O_pe_H(1)];
        elseif abs(O_pe_H(2)-lowest_point(2))<0.02
            coe_planar = [0,1,O_pe_H(2)];
        else
            coe_planar = [O_pe_H(1) O_pe_H(2);lowest_point(1) lowest_point(2)]\[1;1];
            coe_planar(3) = 1;
        end
        
        cy = cos(theta_y_pe_H);
        sy = sin(theta_y_pe_H);
        cx = cos(theta_x_pe_H);
        sx = sin(theta_x_pe_H);
        syms x y
            f11  = (x - O_x_pe_H)*cy-(h_down-O_z_pe_H)*sy;
            f121 = (x - O_x_pe_H)*sx*sy;
            f122 = (y - O_y_pe_H)*cx;
            f123 = (h_down-O_z_pe_H)*sx*cy;
        f1 = f11^2 + (f121 + f122 + f123)^2 - R_peg^2;        
        f2 = coe_planar(1)*x + coe_planar(2)*y - coe_planar(3)*1;
        
        S = solve(f2, f1, 'Real', true);
        
        if isempty(S)
            disp('Equation solve failed!');
            d_lowest = NaN;
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
        plot3(edge_point(1), edge_point(2), edge_point(3), 'mo');hold on
        plot3([edge_point(1) lowest_point(1)],...
              [edge_point(2) lowest_point(2)],...
              [edge_point(3) lowest_point(3)], 'm');
        
        % �ж�����һ���ͶӰ��hole��λ�ù�ϵ
        if abs(edge_point(1)^2 + edge_point(2)^2 - R_hole^2) <= 10^-4 % ����hole�ϣ�����Ӵ�
            disp('three-point contact!');
            flag.contact = 3;
            d_lowest = lowest_point(3) - h_down;
        elseif edge_point(1)^2 + edge_point(2)^2 - R_hole^2 < 0 % ����hole�ڣ�����Ӵ�
            disp('two-point contact!');
            flag.contact = 2;
            d_lowest = lowest_point(3) - h_down;
        elseif edge_point(1)^2 + edge_point(2)^2 - R_hole^2 > 0 % ����hole�⣬�����Ӵ�
            disp('edge-only contact!');
            flag.contact = -1;
            
            d1 = norm(lowest_point(1:2) - edge_point(1:2));
            mid_point = [sum(inter_x)/2  sum(inter_y)/2  sum(inter_z)/2].';
            
            d2_cos = ([0;0] - mid_point(1:2))' * (edge_point(1:2) - mid_point(1:2));
%             d2_temp = norm(mid_point(1:2) - edge_point(1:2)) - norm(mid_point(1:2)) - R_hole;
            if d2_cos >= 0 % mid_point��O_hole��edge_point_proj��ͬ��
                d2 = norm(mid_point(1:2) - edge_point(1:2)) - norm(mid_point(1:2)) - R_hole;
            else % mid_point��O_hole��edge_point_proj���м�
                d2 = norm(mid_point(1:2) - edge_point(1:2)) + norm(mid_point(1:2)) - R_hole;
            end
                        
            d3 = abs(h_down - lowest_point(3));
            d4 = d2*d3/d1;
            
            d_lowest = lowest_point(3) - (h_down - d4);
        end
    else % ����߶����ϴ�˵��ֻ��һ���Ӵ��㣬�˴����ж�һ���ǲ��ǲ���Ӵ�(2014-08-22)        
        px = points(1,:);
        py = points(2,:);
        pz = points(3,:);

        ca = axis_o_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cb = axis_n_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cc = axis_a_pe' * axis_a_pe_H / norm(axis_a_pe_H);
        cs = ca^2 + cb^2;

        cdelta = (px*ca + py*cb).^2 - cs*(px.^2 + py.^2 - R^2);

        idx = cdelta>=0;
        
        px = px(idx);
        py = py(idx);
        pz = pz(idx);
        cdelta = cdelta(idx);

        t1 = (-(px*ca+py*cb) + sqrt(cdelta)) / cs;

        ppx = px + t1 * ca;
        ppy = py + t1 * cb;
        ppz = pz + t1 * cc;

        plot3(ppx, ppy, ppz,'r.');
       
        [~, idx] = min(ppz);
        plot3(ppx(idx), ppy(idx), ppz(idx), 'ko');

        if abs(O_pe_H(1)-lowest_point(1))<10^-6%0.02  % ƽ��ӽ���X���y���غ�ʱ�����̲�̬���������棬��ֱ�Ӽ���
            coe_planar = [1,0,O_pe_H(1)];
        elseif abs(O_pe_H(2)-lowest_point(2))<10^-6%0.02
            coe_planar = [0,1,O_pe_H(2)];
        else
            coe_planar = [O_pe_H(1) O_pe_H(2);lowest_point(1) lowest_point(2)]\[1;1];
            coe_planar(3) = 1;
        end
        
        cd = abs(coe_planar(1)*ppx + coe_planar(2)*ppy - coe_planar(3)*1) < 10^-4;
        [~, idx] = min(ppz(cd));
        
        if isempty(idx) || min(z) < ppz(idx) % һ��Ӵ�            
            disp('one-point contact!');
            flag.contact = 1;

            h_down = min(z);
        else
           disp('edge-only contact!');
            flag.contact = -1;
            
            h_down = ppz(idx);            
        end
        d_lowest = lowest_point(3) - h_down;
    end
end

%% �쳣����
if d_lowest > 5 || d_lowest < -5
    d_lowest = NaN;
    fprintf(2,'Unexpected result!\n');
end

%% ��͵���ֵ
title(['d_{lowest} = ' num2str(d_lowest)]);
disp(flag);





