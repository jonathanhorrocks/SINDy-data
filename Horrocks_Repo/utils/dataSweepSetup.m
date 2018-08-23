function [t, x] = dataSweepSetup(S0, fac, location)

%% Import Data

%Set up initial parameters

alpha = 8.5;
sMethod = 'FGlocal';
iMethod = 'Chris';


switch location
    case 'UK'
    years = 1952:1966;
    pp = 0.95;
    D_i = 2.0;
    L = 65.0*52.0;
    % Imports measles data matrix.
    % Columns: | Time | Cases | Population | Births |
    mData = importdata('mDataEW_N.mat');
    mData = mData(208:end, :); 
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

% Reconstruct Susceptibles

S0 = S0*P(1);
switch sMethod
    case 'constant birth'
        S = SuscRec(C, B, S0, alpha);
        S = S./P;
        
    case 'residuals'
        Sh = sHat(C);
        load enso
        f = fit(t,Sh,'smoothingspline');
        for i = 1:length(t)
            myFit(i) = f(t(i));
        end
        myFit = myFit( : );
        S = abs(myFit - Sh);
        S = S/norm(S, 2);
        
    case 'time delay'
        disp('Run the other script dummy')
        
    case 'FG'
        S(1) = S0;
        alphaV = zeros(N, 1);
        alphaV(1) = alpha;
        Z(1) = 0;
        %hbar = waitbar(0,'Computing remainders...');
        for i = 2:N
            [Z(i), alphaV(i)] = SuscRec_FG(C(1:i), B(1:i), fac);
            S(i) = S0 + Z(i);
            %waitbar(i/N)
        end
        %close(hbar)
        S = S(:);
        S = S./P;
    case 'FGlocal'
        S(1) = S0;
        alphaV = zeros(N, 1);
        alphaV(1) = alpha;
        Z(1) = 0;
        %hbar = waitbar(0,'Computing remainders...');
        for i = 2:N
            [Z(i), alphaV(i)] = SuscRec_FGlocal(C(1:i), B(1:i), fac);
            S(i) = S0 + Z(i);
            %waitbar(i/N)
        end
        %close(hbar)
        S = S(:);
        S = S./P;
    otherwise
        warning('No susceptible reconstruction method selected')
end

% Convert Incidence to Prevalence
switch iMethod
    case 'Simple'
        I = alpha*C;
        I = I./P; 
    case 'Mine'
        I = zeros(N,1);
        if strcmp(location, 'UK')
            I(1) = 3500;
            muI = 1.57;
            sigmaI = .28;
        elseif strcmp(location, 'Ontario')
            I(1) = 400;
            muI = 1.71;
            sigmaI = .14;
        else
            warning('You location-goofed');
        end

        for i = 3:N
            Di = normrnd(muI, sigmaI);
            if Di < 1
                Di = 1;
            end
            iTemp = (1 - (1/Di))*I(i-1) + alpha*C(i);
            if iTemp>0
                I(i) = iTemp;
            iTemp = 0;   
            end           
        end
        I = I./P; 
    case 'Chris'
        cbar = mean(C);
        I = C*pp*D_i/(cbar*L);
end

%Make the data nice and smooooooth (these can and should be changed)       
x(:, 1) = smooth(S, 'sgolay', 2);
x(:, 2) = sgolayfilt(I,3,19);
    