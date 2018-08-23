function [Xi, xA, AIC, SSE] = dataSweepRun(t, x, S0, lambdas, pWs, location, polyorder, folder_path)
usesine = 0; 
useConstant = 1; 
useSeasonal = 1; 
useBirth = 0; %IF YOU TURN THIS ON YOU NEED TO FIX THE ARGUMENTS OF sparsifyDynamics!!
usePSD = 1; %Evaluates model using SSE PSD instead of AIC
omega = 52; %Time-step (bad notation)
n = 2; %Number of state variables
T = 1; %Seasonal period
t_as = 1; %Starting point for analysed time
forcing = 'Seasonal';
dispFigs = 1;

%% Shift Time-Series
x_1 = x(1:end-1, :); %Data matrix
x_2 = x(2:end, :);  %Response matrix

%Start evaluation at later time
x_1 = x_1(t_as:end, :);
x_2 = x_2(t_as:end, :);
Nnew = length(x_1) + 1;
tNew = t(t_as:end, :);

%% Choose Transmission Forcing
switch forcing
    case 'Seasonal'
        peakW = 11; %Timestep of peak transmission 
        betaT = cos(2*pi*(tNew - peakW/52)/T);
        
    case 'Term-Time'
        betaT = zeros(N, 1);
        beta_Yearly = termTimeB(xSpl, P, nY);
        for i = 1:N
            rem = mod(t(i), 1);
            rem_week = int16(52*rem);
            if rem_week == 0
                rem_week = 52;
            end
            betaT(i) = beta_Yearly(rem_week );
        end  
        betaT = betaT/mean(betaT);
end    

%% pool Data  (i.e., build library of nonlinear time series)
if ~useBirth
    Theta = poolDataSeasonal(tNew(1:end-1), x_1,n,polyorder,usesine, 1, omega, betaT(1:end-1));
else
    Theta = poolDataSeasonal(tNew(1:end-1), x_1,n,polyorder,usesine, 1, omega, betaT(1:end-1), B(1:end-1));
end


%% compute Sparse regression: sequential least squares
method = 'GRID';

switch method
    case 'Basic'
        lambda = lambdas;
        Xi(:, 1) = sparsifyDynamics(Theta,x_2(:, 1),lambda, 1);
        Xi(:, 2) = sparsifyDynamics(Theta,x_2(:, 2),lambda, 1);
    case 'AIC'
        deets = [polyorder, usesine, useConstant, useSeasonal, useBirth, T, peakW];
        [Xi, ~, ~, ~, ~] = sparsifyDynamicsAIC_D(Theta,x_2,lambdas,n, deets,x,betaT, 0., t);
    case 'WL'
        allXis = WeakestLink(Theta, x_2, n);
        Xi = allXis(:, 21:22);
    case 'GRID'
        peakWs = pWs;
        bestXis_list = [];
        allXis_list = [];
        allAICs_list = [];
        count = 0;
        for pW = peakWs
            count = count + 1;
            betaT = cos(2*pi*(tNew - pW/52)/T);
            if ~useBirth
                Theta = poolDataSeasonal(tNew(1:end-1), x_1,n,polyorder,usesine, 1, omega, betaT(1:end-1));
            else
                Theta = poolDataSeasonal(tNew(1:end-1), x_1,n,polyorder,usesine, 1, omega, betaT(1:end-1), B(1:end-1));
            end
            deets = [polyorder, usesine, useConstant, useSeasonal, useBirth, T, pW];
            [Xi, ~, allXis, allAICs, ~] = sparsifyDynamicsAIC_Dpsd(Theta,x_2,lambdas,n, deets,x_1,betaT, 0., t);
            bestXis_list = cat(2, bestXis_list, Xi);
            allXis_list = cat(2, allXis_list, allXis);
            allAICs_list = cat(2, allAICs_list, allAICs);
        end
    otherwise
        disp('oops');
end

%% Find min AIC, select that model
if strcmp(method, 'GRID')
    [min_aic, m_ind] = min(allAICs_list);
    Xi = allXis_list(:, m_ind*2 - 1:m_ind*2);
    
    %Get lambda, peakW used for min AIC
    lmb_ind = mod(m_ind, size(lambdas, 2));
    if lmb_ind == 0
        min_lmb = lambdas(end)
    else
        min_lmb = lambdas(lmb_ind)
    end        
    pW_ind = (m_ind - lmb_ind)/size(lambdas, 2);
    if pW_ind == 0
        min_pW = peakWs(end)
    else
        min_pW = peakWs(pW_ind)
    end
end    

%Displays model coefficients
myYout = poolDataLISTsnl({'S','I'},Xi,n,polyorder,usesine, useConstant, useBirth);

%Evaulates model, starting a few years back to account for transient stage
t0 = t(1);
years_back = 0;
trans_stage = years_back*omega;

xA_full(1, :) = x_1(1, :); 
t_start = t0 - years_back;
tA = linspace(t_start, t(end), trans_stage + size(t, 1));
tA = tA(:);
betaTA = cos(2*pi*(tA - min_pW/52)/T);
for i = 2:Nnew + trans_stage
    xA_full(i, :) = modelEvalDiscrete(xA_full(i-1, :), Xi, polyorder, useConstant, useSeasonal, betaTA(i-1, :), useBirth, 0.);
end

xA = xA_full(trans_stage+1:end, :);

%Find AIC and SSE of SINDy model
if usePSD
    AIC = PSD_ls(x(:, 2), xA(:, 2));
else
    AIC = myAIC(x(t_as:end, :), xA, n, x_1(t_as, :), nnz(Xi));
end
SSEv = cumsum((x-xA).^2);
SSE = SSEv(end, 2);

%% Display time-series plots of data and recovered model
q = Nnew;
if dispFigs
    title = strcat(location, '_', 'S0=', num2str(S0), 'lambda=', num2str(min_lmb));
    f1 = figure('Name',strcat(title, '_Susceptible'));
    figure(f1);
    set(f1, 'Visible', 'off');
    plot(t(1:q), x(1:q, 1), t(1:q), xA(1:q, 1));
    Simg_path = strcat(folder_path, '\', location, '_', 'S0=', num2str(S0), 'lambda=', num2str(min_lmb), '_S_img.png');
    saveas(f1, Simg_path);
    f2 = figure('Name', strcat(title, '_Infected'));
    figure(f2)
    set(f2, 'Visible', 'off');
    plot(t(1:q), x(1:q, 2), t(1:q), xA(1:q, 2));
    Iimg_path = strcat(folder_path, '\', location, '_', 'S0=', num2str(S0), 'lambda=', num2str(min_lmb), '_I_img.png');
    saveas(f2, Iimg_path);
end

%% Save data
try
    ts_data = cell(Nnew + 1, 5);
    ts_data(1, :) = [{'Time'}, {'Sd'}, {'Id'}, {'Sa'}, {'Ia'}];
    ts_data(2:Nnew+1, 1) = num2cell(t);
    ts_data(2:Nnew+1, 2) = num2cell(x(:, 1));
    ts_data(2:Nnew+1, 3) = num2cell(x(:, 2));
    ts_data(2:Nnew+1, 4) = num2cell(xA(:, 1));
    ts_data(2:Nnew+1, 5) = num2cell(xA(:, 2));

    data_path = strcat(folder_path, '\', location, '_', 'S0=', num2str(S0), 'lambda=', num2str(min_lmb), '_data.xls');
    model_path = strcat(folder_path, '\', location, '_', 'S0=', num2str(S0), 'lambda=', num2str(min_lmb), '_model.xls');
    xlswrite(data_path,ts_data);
    myYout(1, 1) = {'Term'};
    myYout(1, 2) = {'S Equation'};
    myYout(1, 3) = {'I Equation'};
    index_col = cell(size(myYout, 1), 1);
    index_col(1, 1) = {'Index'};
    index_col(2:end) = num2cell(1:size(myYout, 1)-1);
    myYout = cat(2, myYout, index_col);
catch
    warning('Excel messed up again');
end
xlswrite(model_path, myYout);

