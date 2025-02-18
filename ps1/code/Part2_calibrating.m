%Author: Nan Jiang
%Date: 2/18/2025
%Purpose: Problem Set 1 Part2

clc;
clear;

global theta mu beta lambda cost_c price_c
% load data
load('../data/data_square_world.mat');
E_data    = data.E;           % gross output
X_data    = data.X;           % trade flow 
N         = data.N;
beta      = data.beta;
mu        = data.mu;
lambda    = data.lambda;

% Parameters
theta  =[4,4];                % trade elasticity
sigma  =[3,3];                % elasticity of substituion in two sectors

% the constant term in the input bundle cost equation
cost_c =[beta(1)^(-beta(1))*lambda(1,1)^(-lambda(1,1))*lambda(2,1)^(-lambda(2,1)),beta(2)^(-beta(2))*lambda(1,2)^(-lambda(1,2))*lambda(2,2)^(-lambda(2,2))];

% the constant term in the price index equation
price_c =[gamma(1+(1-sigma(1))/theta(1))^(1/(1-sigma(1))), gamma(1+(1-sigma(2))/theta(2))^(1/(1-sigma(2)))];


%% RECOVER TRADE COSTS
% recover residuals.
I         = length(E_data);                     % recover number of countries
FE_o      = kron(eye(I),ones(I,1));             % make origin FE
FE_d      = kron(ones(I,1),eye(I));             % make destination FE
X_var     = [FE_o, FE_d];                       % EV           

% run gravity 
d         = zeros(size(X_data));
for k =1:2    

  x_k     = X_data(:,:,k);                      % trade flow in sector k
  y_var   = log(x_k(:));                        % DV 
  beta_g    = pinv(X_var'*X_var)*(X_var'*y_var);  % plim to avoid multicollinearity
  res      = y_var - X_var*beta_g;

  % with residuals, we can use an assumption on trade elasticity
  tcost      = exp(res).^(-(1/theta(k)));
  tcost_mat  = reshape(tcost, [I, I]);
  d(:,:,k)   = tcost_mat./diag(tcost_mat);
   
end

%% RECOVER TFPS

% prepare targets in the data
E_target = E_data/E_data(1);

% guesses and tolerances
T               = ones(I,2);  % guess for T
w               = ones(I,1);  % guess for w
adj_exp_N       = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum_N       = 0.5;        % update parameter sum
adj_exp_T       = 0.5;        % update parameter exp (helpful with extreme values)
adj_sum_T       = 0.5;        % update parameter sum

% loop related parameters
dif_tol_outer   =1;
dif_tol_T       = 1;          % initial tolerance
tol_T           = 1e-6;       % tolerance parameter
count_T         = 1;          % count
max_count_T     = 1000;       % max count
update_outer    = 1;
options = optimoptions('lsqnonlin', 'Display', 'off'); % Suppress output

% search for values of T in the outer loop
while dif_tol_outer > tol_T && count_T < max_count_T && update_outer > tol_T

   
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
    P_sol   = lsqnonlin(F_price, P_guess, lb, [], options); 
    
    % inferred expenditure share based on wage and price index
    pi      = exp_share(P_sol,w,T,d);  
    
    % solve the total expenditure for each country in each sector
    F_expenditure  = @(X) total_expenditure(X,pi,E);
    X_sol          = lsqnonlin(F_expenditure, X_guess, lb, [], options);  %gross output
    
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

end

    % compare model implied with actual gross output
    E_model   = X_sol;
    E_model   = E_model./E_model(1);
    dif_adj_T = E_target./E_model;
    dif_tol_T = max(abs(1-dif_adj_T));
    T_new     = T.*dif_adj_T.^(adj_exp_T);
    T_new     = T_new./T_new(1);    
    T         = T*adj_sum_T + T_new*(1-adj_sum_T);
    T         = T/T(1); % ensure normalization
    update_outer  = abs(max(dif_tol_T)-dif_tol_outer);
    dif_tol_outer = max(dif_tol_T);
    
    

    % display results
    disp(['Tolerance in N (inner) = ', num2str(dif_tol)]); 
    disp(['Tolerance in T (outer) = ', num2str(dif_tol_outer)]); 
    disp(['count (outer) = ', num2str(count_T)]);
    count_T   = count_T + 1; 
end

% recover trade flow, first dimension is the origin, the second dimension
% is the destination, third is the sector
X        = zeros(size(pi));
X(:,:,1) = pi(:,:,1) .*X_sol(:,1)';
X(:,:,2) = pi(:,:,2) .*X_sol(:,2)';

calibration = struct(...
        'd', d, ...
        'T', T, ...
        'N', N, ...
        'I', I);
save('../data/calibration.mat', 'calibration');
