function T = trans(vector)
% Perform translation transformation
%
% Parameters:
%  - vector: translation vector
%
% Author:
%  Ray Lee(lirui2013@ia.ac.cn)
%
% Date:
%  2014-06-14 File created
%  2014-07-17 Add comments
%

[m, n] = size(vector);
if ~((m==3&&n==1)||(m==1&&n==3))
    disp('error!');
    return;
end

if (m<n)
    vector = vector';
end

R = eye(3);
T = [R vector;0 0 0 1];