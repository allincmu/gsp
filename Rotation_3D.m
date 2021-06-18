close all; clear all; clc


N = 3;
ROT_X = 0
ROT_Y = 0
ROT_Z = pi/2

% Generate random N x N matrix
% mtx = rand(N,N)
mtx = rand(3,3);
while rank(mtx) ~= N 
  mtx = rand(3,3);
end

for i = 1:3
  mtx = orth(mtx);
end

% mtx = [1,1,1; 0,0,0; 2,2,2]

% Generate rotation matrix
syms t
syms u
syms v

Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = [cos(u) 0 sin(u); 0 1 0; -sin(u) 0 cos(u)];
Rz = [cos(v) -sin(v) 0; sin(v) cos(v) 0; 0 0 1];

R = Rx*Ry*Rz;
R_rot = subs(R, {t,u,v}, {ROT_X, ROT_Y, ROT_Z});

mtx_rotated = single(subs(R, {t,u,v}, {ROT_X, ROT_Y, ROT_Z})*mtx(1));
mtx_rotated = single(subs(R, {t,u,v}, {ROT_X, ROT_Y, ROT_Z})*mtx(2));
mtx_rotated = single(subs(R, {t,u,v}, {ROT_X, ROT_Y, ROT_Z})*mtx)

a11 = acos(dot(mtx_rotated(:,1), mtx(:,1))/ ( norm(mtx_rotated(:,1)) * norm(mtx(:,1))))
a12 = acos(dot(mtx_rotated(:,1), mtx(:,2))/ ( norm(mtx_rotated(:,1)) * norm(mtx(:,2))))
a13 = acos(dot(mtx_rotated(:,1), mtx(:,3))/ ( norm(mtx_rotated(:,1)) * norm(mtx(:,3))))

a22 = acos(dot(mtx_rotated(:,2), mtx(:,2))/ ( norm(mtx_rotated(:,2)) * norm(mtx(:,2))))
a21 = acos(dot(mtx_rotated(:,2), mtx(:,1))/ ( norm(mtx_rotated(:,1)) * norm(mtx(:,1))))
a23 = acos(dot(mtx_rotated(:,2), mtx(:,3))/ ( norm(mtx_rotated(:,3)) * norm(mtx(:,3))))

a33 = acos(dot(mtx_rotated(:,3), mtx(:,3))/ ( norm(mtx_rotated(:,3)) * norm(mtx(:,3))))
a31 = acos(dot(mtx_rotated(:,3), mtx(:,1))/ ( norm(mtx_rotated(:,3)) * norm(mtx(:,1))))
a32 = acos(dot(mtx_rotated(:,3), mtx(:,2))/ ( norm(mtx_rotated(:,3)) * norm(mtx(:,2))))

figure;
origin = [0,0,0];
p = [origin;mtx(:,1)';nan(1,3);origin;mtx(:,2)';nan(1,3);origin;mtx(:,3)']
plot3(p(:,1),p(:,2),p(:,3))
hold on
p = [origin;mtx_rotated(:,1)';nan(1,3);origin;mtx_rotated(:,2)';nan(1,3);origin;mtx_rotated(:,3)']
plot3(p(:,1),p(:,2),p(:,3))
grid on


% plotv(mtx_rotated, 'r')

% figure
% [px, py, pz] = find_p(mtx);
% surf(px, py, pz)
% 
% figure
% px, py, pz = find_p(mtx_rotated);
% surf(px, py, pz)
% 
% function [px, py, pz] =  find_p(mtx)
% 
%   v1 = mtx(:,1)';
%   v2 = mtx(:,2)';
%   v3 = mtx(:,3)';
%   
%   p1=@(q,r,s) v1(1)*q+v2(1)*r+v3(1)*s;
%   p2=@(q,r,s) v1(2)*q+v2(2)*r+v3(1)*s;
%   p3=@(q,r,s) v1(3)*q+v2(3)*r+v3(1)*s;
%   
%   [Q, R, S]=meshgrid(-5:1:5);
%    
%   px=p1(Q, R, S);
%   py=p2(Q, R, S);
%   pz=p3(Q, R, S);
% 
% end
