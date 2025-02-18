function pi=exp_share(price,wage,T,d)

global theta beta lambda cost_c 

cost = zeros(size(price));        % cost of an input bundle

 for i=1:100
     for k=1:2
         cost(i,k)= cost_c(k)*wage(i)^beta(k)*price(i,1)^lambda(1,k)*price(i,2)^lambda(2,k);
     end
 end

 Pi = zeros(size(price));         % Multilateral Risistance term
 pi = zeros(size(d));             % expenditure share
 % expenditure shares
 for k=1:2
    pi_nom  = T(:,k).*(d(:,:,k).*cost(:,k)).^(-theta(k)); % dim 1 i origin - dim 2 is destination
    Pi(:,k)      = sum(pi_nom,1);       % sum over dim 1 
    pi(:,:,k)      = pi_nom./Pi(:,k);
 end