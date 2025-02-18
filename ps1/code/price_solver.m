%This function caclaue the price index based on wage guess
function F=price_solver(wage,price,T,d)


global theta beta lambda cost_c price_c

cost = zeros(size(price));        % cost of an input bundle

 for i=1:100
     for k=1:2
         cost(i,k)= cost_c(k)*wage(i)^beta(k)*price(i,1)^lambda(1,k)*price(i,2)^lambda(2,k);
     end
 end

 Pi = zeros(size(price));         % Multilateral Risistance term
 % expenditure shares
 for k=1:2
    pi_nom  = T(:,k).*(d(:,:,k).*cost(:,k)).^(-theta(k)); % dim 1 i origin - dim 2 is destination
    Pi(:,k)      = sum(pi_nom,1);       % sum over dim 1
 end

F=zeros(200,1);
for i=1:100
    for k =1:2
        F((i-1)*2+k,1)=price(i,k)-price_c(k)*Pi(i,k)^(-1/theta(k));
    end
end
end