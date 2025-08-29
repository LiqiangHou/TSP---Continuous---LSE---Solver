
%%
function [J ,Y] = path_cost_X(N, X)

dimension = 2;
for i=1:N
    X_i = X((i - 1)*dimension+1 : (i - 1)*dimension+2);

    Y(i) = norm(X_i);
end
J = sum(Y);
return
end


