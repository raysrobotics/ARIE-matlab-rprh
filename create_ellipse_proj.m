function P = create_ellipse_proj(r, theta, transform)
% ͨ����Բ���̽�����ת�õ�����淽��
%
%   r - radius of the circle
%   theta - parameter of the parameter equation of the circle
%           x = r*cos(theta)
%           y = r*sin(theta)
%   transform - transformation matrix T
%
xx = r*cos(theta);
yy = r*sin(theta);
zz = 0*theta;
ww = 0*theta+1;

P = transform * [xx;yy;zz;ww];
