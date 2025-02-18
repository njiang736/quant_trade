%this function calculate the total expenditure 
function F=total_expenditure(X,pi,E)

global mu lambda

F=zeros(200,1);
for i=1:100
    for k=1:2
        F((i-1)*2+k,1)= X(i,k)-lambda(k,1)*pi(i,:,1)*X(:,1)-lambda(k,2)*pi(i,:,2)*X(:,2)-mu(k)*E(i);
    end
end
end