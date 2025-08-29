
function [bar_X,W,R,M_T] = Init_X_Set(graph)
global R
global T
%
init_step = [ graph.min_abs_dst_x  , ...
              graph.min_abs_dst_y  ];

min_d_pos     = [graph.min_dx        graph.min_dy];
max_d_pos     = [graph.max_dx        graph.max_dy];
max_std_dr    = [graph.max_abs_dst_x  ,...
    graph.max_abs_dst_y  ];
%

bar_X = [];
for i=1:graph.n 
    % bar_X = [bar_X , -5.0*init_step];
    bar_X = [bar_X , 1.0*init_step];
end
bar_X = bar_X';
%

sigma_x0 = [1.0, 1.0];

init_step = [ graph.min_abs_dst_x /(graph.max_dx - graph.min_dx) , ...
    graph.min_abs_dst_y /(graph.max_dy - graph.min_dy) ] ;
mu_x0     = 0.5 - init_step;
 

% 



for i=1:graph.n
    Q{i} =  [  graph.min_abs_dst_x   ,  0.0
               0.0                   ,  graph.min_abs_dst_y ] *150;     % 
end



% T= diag([1.0 1.0].^2);
T= diag([graph.min_abs_dst_x graph.min_abs_dst_y].^2);


% 
for i=1:graph.n
    W(((i-1)*2 + 1):((i-1)*2 + 2) , ((i-1)*2 + 1):((i-1)*2 + 2)) = Q{i};  
    M_T(((i-1)*2 + 1):((i-1)*2 + 2) , ((i-1)*2 + 1):((i-1)*2 + 2)) = T; 
end

% 
R = (graph.min_edge)^2;              
% 


return
end

