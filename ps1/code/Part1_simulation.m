%Author: Nan Jiang
%Date: 2/14/2025
%Purpose: Problem Set 1 Part1

clc;
clear;

% load fundamentals
load("..\data\fundamental_for_simulation.mat")
d     = fundamentals.d; 
I     = fundamentals.I;
T     = fundamentals.T;
N     = fundamentals.N;
coord = fundamentals.coord;

% Parameters
global theta mu beta lambda cost_c price_c
theta  =[4,4];                % trade elasticity
mu     =[0.2,0.8];            % consumption share in each sector
beta   =[0.5,0.5];            % cost share of labor
lambda =[0.4,0.1;0.1,0.4];    % cost share of composite goods
sigma  =[3,3];                % elasticity of substituion in two sectors

% the constant term in the input bundle cost equation
cost_c =[beta(1)^(-beta(1))*lambda(1,1)^(-lambda(1,1))*lambda(2,1)^(-lambda(2,1)),beta(2)^(-beta(2))*lambda(1,2)^(-lambda(1,2))*lambda(2,2)^(-lambda(2,2))];

% the constant term in the price index equation
price_c =[gamma(1+(1-sigma(1))/theta(1))^(1/(1-sigma(1))), gamma(1+(1-sigma(2))/theta(2))^(1/(1-sigma(2)))];

% solve the model
w    = ones(I,1);             % initial guess for wage
dif_tol   = 1;                % initial tolerance
tol       = 1e-6;             % tolerance parameter
count     = 1;                % count
max_count = 1000;             % max count
update    = 1;                % difference of the distance in erros bewteen each loop
adj_exp   = 0.5;              % update parameter exp (helpful with extreme values)
adj_sum   = 0.5;              % update parameter sum
P_guess   = ones(I,2);                       % initial guess for price index
lb        = zeros(size(P_guess));            % non negative constraint
X_guess   = ones(I,2);                       % initial guess for total expenditure


while dif_tol > tol && count < max_count && update > tol

    % income
    w       = w./w(1);  % ensure first country is numeraire  
    E       = w.*N;
    
    % solve the price index based on the wage guess
    F_price = @(price) price_solver(w,price,T,d);
    P_sol   = lsqnonlin(F_price, P_guess, lb, []); 
    
    % inferred expenditure share based on wage and price index
    pi      = exp_share(P_sol,w,T,d);  
    
    % solve the total expenditure for each country in each sector
    F_expenditure  = @(X) total_expenditure(X,pi,E);
    X_sol          = lsqnonlin(F_expenditure, X_guess, lb, []);
    
    % demand for labor
    N_dem     = N_demand(X_sol,pi,w);
    dif_adj_N = N_dem./N;
    dif_tol_N = max(abs(1-dif_adj_N));
    w_new     = w.*dif_adj_N.^adj_exp;
    w_new     = w_new./w_new(1); 
    w         = w*adj_sum + w_new*(1-adj_sum); 

    % update
    update  = abs(dif_tol_N-dif_tol);
    dif_tol = dif_tol_N;
    count   = count + 1;

    % display results
    disp(['Tolerance in N = ', num2str(dif_tol_N)]) 
    disp(['Iteration = ', num2str(count)]) 


end

% recover trade flow, first dimension is the origin, the second dimension
% is the destination, third is the sector
X        = zeros(size(pi));
X(:,:,1) = pi(:,:,1) .*X_sol(:,1)';
X(:,:,2) = pi(:,:,2) .*X_sol(:,2)';

% export data
data = struct(...
        'N',N, ...
        'E',X_sol, ...  % gross output
        'X', X, ...     % tradee flow
        'mu', mu, ...
        'beta', beta,...
        'lambda',lambda);
save('../data/data_square_world.mat', 'data');
