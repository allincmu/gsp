
M = rand(10);
M = (M + M') ./ 2;

[V, D] = eig(M);

for i = 1:length(V)
  M*V(:,i) - D(i)*V(:,i)
end