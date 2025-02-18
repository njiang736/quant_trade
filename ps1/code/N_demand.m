function N=N_demand(X,pi,w)

global beta

N=zeros(size(w));

for i=1:100

    N(i) = (beta(1)*pi(i,:,1)*X(:,1)+ beta(2)*pi(i,:,2)*X(:,2))/w(i);

end