%%

function [F] = path_Jacobian(N, X)

F   = [];        % The Jacobian in vector form
dimension = 2;


for i = 1:N
    x_i = X((i - 1)*dimension + 1 : (i - 1)*dimension + 2);
    f_i = [x_i(1)/sqrt(x_i(1)^2 + x_i(2)^2)  , x_i(2)/sqrt(x_i(1)^2 + x_i(2)^2)];

    F   = [F , f_i];

end

return
end