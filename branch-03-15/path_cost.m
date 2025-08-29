%%

function [J , Y] = path_cost(I_path, edge)

for i=1:length(I_path)-1
    Y(i)= edge(I_path(i) , I_path(i+1)); 
end

J = sum(Y);
return
end
