function [sigma] = SourceStrengths(n_hat,U_infty,alpha)

Ni = length(n_hat(1,:));
sigma = zeros(Ni,1);

for i = 1:Ni
    
    sigma(i) = U_infty*[cos(alpha),sin(alpha)]*n_hat(:,i);
    
end