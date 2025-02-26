%Author: Nan Jiang
%Date: 2/18/2025
%Purpose: Problem Set 1 Part3

clc;
clear;

% load data
load('../data/calibration.mat', 'calibration');
d     = calibration.d;
T     = calibration.T;
N     = calibration.N;
I     = calibration.I;

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

%% Counterfactual
w    = ones(I,1);             % initial guess for wage
adj_exp   = 0.5;              % update parameter exp (helpful with extreme values)
adj_sum   = 0.5;              % update parameter sum
P_guess   = ones(I,2);                       % initial guess for price index
lb        = zeros(size(P_guess));            % non negative constraint
X_guess   = ones(I,2);                       % initial guess for total expenditure
options = optimoptions('lsqnonlin', 'Display', 'off'); % Suppress output

% we will run this for three Ts. We will increase T(1,1) by 20 percent for the
% second round, and increase T(2,1) by 20 percent for the thrid round
T_case      = [T,T,T];
T_case(1,3) = T_case(1,3)*1.2;
T_case(2,5) = T_case(2,5)*1.2;
W_case      = zeros(I,3);                    % welfare in each time
w_case      = zeros(I,3);                    % wage in each time

for tx =1:3

    % pick T    
    T = T_case(:,tx*2-1:tx*2);

    dif_tol   = 1;                % initial tolerance
    tol       = 1e-6;             % tolerance parameter
    count     = 1;                % count
    max_count = 1000;             % max count
    update    = 1;                % difference of the distance in erros bewteen each loop

    while dif_tol > tol && count < max_count && update > tol

    % income
    w       = w./w(1);  % ensure first country is numeraire  
    E       = w.*N;
    
    % solve the price index based on the wage guess
    F_price = @(price) price_solver(w,price,T,d);
    P_sol   = lsqnonlin(F_price, P_guess, lb, [],options); 
    
    % inferred expenditure share based on wage and price index
    pi      = exp_share(P_sol,w,T,d);  
    
    % solve the total expenditure for each country in each sector
    F_expenditure  = @(X) total_expenditure(X,pi,E);
    X_sol          = lsqnonlin(F_expenditure, X_guess, lb, [],options);
    
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

    % display results
    disp(['Tolerance in N = ', num2str(dif_tol_N)]) 
    disp(['Iteration = ', num2str(count)]) 
    count   = count + 1;

end
   w_case(:,tx) = w;
   P_index      = (P_sol(:,1)./mu(1)).^mu(1) .* (P_sol(:,2)./mu(2)).^mu(2);
   W_case(:,tx) = w./P_index;

end

%% make some figures

% increase tfp of country 1 in sector 1
D_rincome = 100*(W_case(:,2)-W_case(:,1))./W_case(:,1);
D_rincome(1) = max(D_rincome(2:end)); % FOR VISUALIZATION WE WILL REMOVE THE GAINS FOR 1
I_side  = sqrt(I);
D_rincome_mat = reshape(D_rincome, [I_side, I_side]);
imagesc(D_rincome_mat);
colorbar
saveas(gcf, '../output/map_D_real_income_levels_case1.jpg');

% increase tfp of country 2 in sector 1
D_rincome = 100*(W_case(:,3)-W_case(:,1))./W_case(:,1);
D_rincome(2) = max(D_rincome(3:end)); % FOR VISUALIZATION WE WILL REMOVE THE GAINS FOR 1
I_side  = sqrt(I);
D_rincome_mat = reshape(D_rincome, [I_side, I_side]);
imagesc(D_rincome_mat);
colorbar
saveas(gcf, '../output/map_D_real_income_levels_case2.jpg');

