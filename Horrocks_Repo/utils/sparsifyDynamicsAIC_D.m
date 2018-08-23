%Runs a sparse linear regression solving for the coefficients of
%x_{t+1} = f(x+t)

function [minXi, minAIC, Models, IC, varAIC] = sparsifyDynamicsAIC_D(Theta,dXdt,lambdas,n, deets, x, betaT, Bir, t)
if nargin < 9 || isempty(t)
    t = linspace(0, 0);
end
if nargin < 8 || isempty( Bir )
    Bir = 0.;
end
polyorder = deets(1);
useConstant = deets(3);
useSeasonal = deets(4);
useBirth = deets(5);
T = deets(6);
peakW = deets(7);
useWeighted = 1;
usePresim = 0;
tol = -5.5e17;
t0 = t(1);
t_end = t(end);
years_back = 0;
trans_stage = years_back*T;
S_constraint = 1;

Models = [];
IC = [];
varAIC = [];
for lambda = lambdas    
    % compute Sparse regression: sequential least squares
    Xi = Theta\dXdt;  % initial guess: Least-squares

    % lambda is our sparsification knob.
    if ~useWeighted
        count = 1;
        A = 1000;
        
        %Continues sparsifying until AIC is low enough (or count gets too
        %high)
        max = 10;
        hbar = waitbar(0,'Making those dynamics extra sparse...');
        while A > tol && count < max
            waitbar(count/max)
            count = count + 1;
            smallinds = (abs(Xi)<lambda);   % find small coefficients
            Xi(smallinds)=0;                % and threshold
            
            %Regresses dynamics for each state variable
            for ind = 1:n                   % n is state dimension
                biginds = ~smallinds(:,ind);
                % Regress dynamics onto remaining terms to find sparse Xi
                Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind); 
            end
        end
        close(hbar)
    else %Uses weighted lambda, without AIC stopping
        for j = 1:10
            smallinds = zeros(size(Xi, 1), n);
            for k = 1:size(Xi, 1)
                w = norm(Theta(:, k));
                for v = 1:n
                    if abs(Xi(k, v)) < lambda/w;
                        smallinds(k, v) = 1;
                    end
                end
            end
            smallinds = logical(smallinds);
            Xi(smallinds)=0;                % and threshold
            
            if ~S_constraint
                for ind = 1:n                   % n is state dimension
                    biginds = ~smallinds(:,ind);
                    % Regress dynamics onto remaining terms to find sparse Xi
                    Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind); 
                end
            else %Force coefficients to be opposite and equal, switch each time
                if rem(j, 2)
                    ind = 2;
                    oth = 1;
                else
                    ind = 1;
                    oth = 2;
                end 
                
                % Find lower and upper bound (equal to negative of other
                % vector)
                biginds = ~smallinds(:, ind);
                Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
                
                options = optimoptions('lsqlin','Algorithm','interior-point', 'Display', 'off');
                lb = zeros(nnz(biginds), 1);
                ub = zeros(nnz(biginds), 1);
                cnt = 0;
                if polyorder == 2
                    fo_lb = 6;
                    fo_up = 10;
                elseif polyorder == 3
                    fo_lb = 10;
                    fo_up = 14;
                else
                    fo_lb = 10;
                    fo_up = 14;
                end
                for q = 1:size(biginds, 1)                   
                    if biginds(q) == 1
                        cnt = cnt +1;
                        if q < 4 || ((fo_lb < q) && (q < fo_up))                       
                            lb(cnt) = -Inf;
                            ub(cnt) = Inf;
                        else
                            lb(cnt) = -Xi(q, ind);
                            ub(cnt) = -Xi(q, ind);
                        end
                    end
                end
                
                Xi(smallinds(:, ind), oth) = 0.;
                try
                    Xi(biginds, oth) = lsqlin(Theta(:, biginds), dXdt(:, oth), [], [], [], [], lb, ub, [], options);
                catch e
                    fprintf(1,e.message);
                end
%                 Xi(biginds,oth) = Theta(:,biginds)\dXdt(:,oth);
%                 Xi(biginds(4:end), oth) = -Xi(biginds(4:end), ind);
            end
        end
    end
    % Simulate model, compute AIC
    xA(1, :) = x(1, :);
    if usePresim
        t_start = t0 - years_back;
        tA = linspace(t_start, t_end, trans_stage + size(t, 1));
        tA = tA(:);
        betaT = cos(2*pi*(tA - peakW/52)/T);
    else
        trans_stage = 0;
    end
    for i = 2:length(x) + trans_stage;
        xA(i, :) = modelEvalDiscrete(xA(i-1, :), Xi, polyorder, useConstant, useSeasonal, betaT(i-1), useBirth, Bir);
    end
    xA = xA(trans_stage+1:end, :);
    k = nnz(Xi);
    Models = cat(2, Models, Xi);
    if size(x, 1) == size(xA, 1)
        A = myAIC(x, xA, n, x(1, :), k);
    else
        A = 1000000; %Large number, indicates something went wrong
    end
    IC = cat(2, IC, A);
end

%Sort the AICs, return the model at the smallest one.  
[vals, inds] = sort(IC, 2);
minXi = Models(:, n*(inds(1)-1)+1:n*inds(1));
minAIC = vals(1);