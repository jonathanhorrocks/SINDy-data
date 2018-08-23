function [Zt, alpha] = SuscRec_FGlocal(C, B, fac)
%Returns the weighted difference from susceptible mean
%(Zt-Z0) and the estimated under-reporting rate

Y = cumsum(B);
X = cumsum(C);
rHat = gaussKE(fac*std(X), X ,Y);


alpha = rHat(end)/X(end);
Zt = Y(end) -  rHat(end);