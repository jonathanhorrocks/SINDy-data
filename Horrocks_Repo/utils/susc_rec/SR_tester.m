%Tests Susceptible Reconstruction method.  
location = 'Ontario_Rubella';
f = 0.45;
S0_p = 0.15;

switch location
    case 'UK'
    years = 1948:1966;
    pp = 0.95;
    D_i = 2.0;
    L = 65.0*52.0;
    % Imports measles data matrix.
    % Columns: | Time | Cases | Population | Births |
    mData = importdata('mDataEW_N.mat');
    [row, col] = size(mData);
    t = mData(:, 1);
    C = mData(:, 2);
    P = mData(:, 3);
    B = mData(:, 4);
    N = length(t);
    dt = 1/N;
        
    case 'Ontario'
    years = 1946:1966;
    pp = 0.95;
    D_i = 2.;
    L = 65*52.;
    
    filename = 'OntarioChickenWeekly39_69.txt';
    % delimiterIn = ' ';
    % headerlinesIn = 1;
    cData = importdata(filename);
    cData = cData(366:end, :);
    t = cData(:, 1);
    C = cData(:, 2);
    %C = smooth(C);
    N = length(t);
    dt = 1/N;

    tempData = importdata('Ontario_Demographics_Measles.txt');
    tTemp = tempData(:, 1);
    pTemp = tempData(:, 2);
    P = interp1(tTemp, pTemp, t);

    %Import Birth Rate data
    B_Data = importdata('Ontario_Birth_Data_M.txt');
    tYears = B_Data(:, 1);
    Btemp = B_Data(:, 2);
    B = interp1(tYears, Btemp, t)/(52/4);
    BP = B./P;
    
    case 'Ontario_Rubella'
    years = 1946:1960;
    pp = 0.95;
    D_i = 2.;
    L = 65*52.;
    
    filename = 'OntarioRubellaWeekly39_69.txt';
    % delimiterIn = ' ';
    % headerlinesIn = 1;
    cData = importdata(filename);
    cData = cData(366:1135, :);
    t = cData(:, 1);
    C = cData(:, 2);
    C(33) = 5;
    C(348) = 7;
    %C = smooth(C);
    N = length(t);
    dt = 1/N;

    tempData = importdata('Ontario_Demographics_Measles.txt');
    tTemp = tempData(:, 1);
    pTemp = tempData(:, 2);
    P = interp1(tTemp, pTemp, t);

    %Import Birth Rate data
    B_Data = importdata('Ontario_Birth_Data_M.txt');
    tYears = B_Data(:, 1);
    Btemp = B_Data(:, 2);
    B = interp1(tYears, Btemp, t)/(52/4);
    BP = B./P;
    
    otherwise
        warning('No location selected')
end
    
% 
% filename = 'OntarioChickenWeekly39_69.txt';
% cData = importdata(filename);
% cData = cData(366:end, :);
% t = cData(:, 1);
% C = cData(:, 2);
% %C = smooth(C);
% 
% tempData = importdata('Ontario_Demographics_Measles.txt');
% tTemp = tempData(:, 1);
% pTemp = tempData(:, 2);
% P = interp1(tTemp, pTemp, t);
% 
% %Import Birth data
% B_Data = importdata('Ontario_Birth_Data_M.txt');
% tYears = B_Data(:, 1);
% Btemp = B_Data(:, 2);
% B = interp1(tYears, Btemp, t)/(52/4);

%Reconstruct susceptible class
N = length(t);
S0 = S0_p*P(1);

S(1) = S0;
alpha(1) = 8;
Z(1) = 0;
hbar = waitbar(0,'Computing remainders...');
for i = 2:N
    %[Z(i), alpha(i)] = SuscRec_FG(C(1:i), B(1:i));
    [Z(i), alpha(i)] = SuscRec_FGlocal(C(1:i), B(1:i), f);
    S(i) = S0 + Z(i);
    waitbar(i/N)
end
close(hbar)
S = S(:);
S = S./P;


X = cumsum(C);
Y = cumsum(B);


%  [h, S1, S2] = bwEstimator(C, B); 
%  plot(h, S1, h, S2)

% rHat = gaussKE(0.5*std(X), X, Y);
% P = polyfit(X,Y,1);
% Yhat = P(1)*X + P(2);
%scatter(X, Y);

% for i = 1:N
%     Xvec = cumsum(C(1:i));
%     X = Xvec(end);
%     Yvec = cumsum(B(1:i));
%     Y = Yvec(end);
%     f = @(x) alpha(i)*x + Z(i);
%     scatter(Xvec, Yvec)
% %     hold on;
% %     ezplot(f, Xvec(1), Xvec(i))
% end
    
plot(t, S)