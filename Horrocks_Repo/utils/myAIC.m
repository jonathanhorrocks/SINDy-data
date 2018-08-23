function v = myAIC(x, xM, n, x0, k)

m = length(x);
E = zeros(m, 1);

for i = 1:m
    for j = 1:n
        E(i) = E(i) + abs(x(i, j) - xM(i, j));
    end
end

Esum = cumsum(E);
EsumF = Esum(m);
v = m*log(EsumF/m) + 2*k;

%Correct for finite sample size
v = v + 2*(k+1)*(k+2)/(m - k - 2);