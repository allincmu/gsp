

% Generate random N x N matrix
% mtx = rand(N,N)
mtx = rand(2,2);
while rank(mtx) ~= 2
  mtx = rand(2,2);
end


mtx = orth(mtx);


theta = pi/3;
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

R_rot = R*mtx;

acos(dot(R_rot(:,1), mtx(:,1)))
acos(dot(R_rot(:,2), mtx(:,2)))

