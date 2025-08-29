
%%

function L_i = Likelihood_of_node(N,X,bar_X_i,SIGMA_i_i)
global T

% system noise
% T       = diag([1.0  , 1.0].^2);

for j=1:N

    L_i(j) = (X(j,:) - bar_X_i')*pinv(SIGMA_i_i + T+ eye(size(SIGMA_i_i))*1.0e-10 , 1.0e-20) *(X(j,:) - bar_X_i')' + log(det(SIGMA_i_i + T + eye(size(SIGMA_i_i))*1.0e-10));    


end

return
end


