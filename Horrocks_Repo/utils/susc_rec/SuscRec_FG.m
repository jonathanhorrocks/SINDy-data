function [Zt, alpha] = SuscRec_FG(C, B)
%Returns the weighted difference from susceptible mean
%(Zt-Z0) and the estimated under-reporting rate

Y = cumsum(B);
X = cumsum(C);
P = polyfit(X,Y,1);

alpha = P(1);
Zt = Y(end) -  alpha*X(end);


