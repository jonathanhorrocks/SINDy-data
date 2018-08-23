%Runs the SINDy algorithm over a parameter plane of S_0, phi, and lambda
%values

lambdas = logspace(-4, 0, 15);
pWs = linspace(0, 51.5, 104);
S0s = linspace(0.05, 0.15, 15);
facs = linspace(0.15, 0.45, 15);
%S0s = 0.12;
folder_path = 'output path goes here';
location = 'Ontario_Rubella'; %UK, Ontario, or Ontario_Rubella
polyord = 2;

Xis = [];
xAs = [];
AICs = zeros(size(lambdas, 2), size(S0s, 2));
SSEs = zeros(size(lambdas, 2), size(S0s, 2));
classes = zeros(size(lambdas, 2), size(S0s, 2));
hbar = waitbar(0, 'Running sweep...'); 
for i = 1:size(S0s, 2) 
    S0 = S0s(i);
    fac = facs(i);
    [t, x] = dataSweepSetup(S0, fac, location);
    for j = 1:size(lambdas, 2)
        waitbar(((i-1)*size(lambdas, 2) + j)/(size(S0s, 2)*size(lambdas, 2)));
        lambda = lambdas(j);
        [Xi, xA, AIC, SSE] = dataSweepRun(t, x, S0, lambda, pWs, location, polyord, folder_path);
        
        % Get AIC and SSE values, add vectors to list
        AICs(j, i) = AIC;
        SSEs(j, i) = SSE;
        Xis = cat(2, Xis, Xi);
        xAs = cat(2, xAs, xA);
        
        % Get class (period) of ts, add to list
        cl = getpeaks(xA(:, 2));
        classes(j, i) = cl;
    end
end
close(hbar);

AIC_data = cell(size(AICs)+1);
AIC_data(1, 1) = {'Lambda/S0'};
AIC_data(2:end, 1) = num2cell(lambdas);
AIC_data(1, 2:end) = num2cell(S0s);
AIC_data(2:end, 2:end) = num2cell(AICs);
AIC_path = strcat(folder_path, '\', location, '_PSDs.xls');
xlswrite(AIC_path, AIC_data)

SSE_data = cell(size(SSEs)+1);
SSE_data(1, 1) = {'Lambda/S0'};
SSE_data(2:end, 1) = num2cell(lambdas);
SSE_data(1, 2:end) = num2cell(S0s);
SSE_data(2:end, 2:end) = num2cell(SSEs);
SSE_path = strcat(folder_path, '\', location, '_SSEs.xls');
xlswrite(SSE_path, SSE_data)

CLS_data = cell(size(classes)+1);
CLS_data(1, 1) = {'Lambda/S0'};
CLS_data(2:end, 1) = num2cell(lambdas);
CLS_data(1, 2:end) = num2cell(S0s);
CLS_data(2:end, 2:end) = num2cell(classes);
CLS_path = strcat(folder_path, '\', location, '_classes.xls');
xlswrite(CLS_path, CLS_data)