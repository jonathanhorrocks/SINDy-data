function S = SuscRec(C, B, S0, alpha)

N = length(C);
S = zeros([1, N]);
S(1) = S0;
for i = 2:length(C)
    Stemp = S(i-1) + B(i) - alpha*C(i);
    if Stemp > 0
        S(i) = Stemp;
    else
        S(i) = 0;
    end
end
S = S( : );