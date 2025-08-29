

function [J,X_path,V,I_path] = cost_fun_0517_2(bar_X,F_path,W ,R, graph)
global Q_iminus
global bar_J
global Q_add
global Q_cov
global M_V
global F_likelihood
global F_likelihood_mu_Y



edge = graph.edges;
N    = graph.n;


for i=1:N
    pre_Q{i} = W(((i-1)*2 + 1):((i-1)*2 + 2) , ((i-1)*2 + 1):((i-1)*2 + 2));     % 
end

Q_cov{1} = zeros(2,2);
for i=1:N-1
    Q_cov{i} = W((i*2 + 1):(i*2 + 2) , ((i-1)*2 + 1):((i-1)*2 + 2));     % 
end



% [F_bar]   = path_Jacobian(N, bar_X);

[V,Q] = path_prior_cov(N, F_path,pre_Q,R);



[X_path,I_path] = select(N,bar_X,Q,graph);

%
[bar_J , bar_Y] = path_cost_X(N, bar_X);
[J , Y]         = path_cost(I_path, edge);


for i=1:N
    F_i            = F_path((i - 1)*2+1 : (i-1)*2+2);



    F_likelihood_i = (1 + (Y(i) - bar_Y(i))/M_V(i))*F_i;
    F_likelihood_mu_Y_i = (1 + (Y(i) - bar_Y(i))/(M_V(i) + bar_Y(i)^2))*F_i;          % varaiance oabtined using convoultion 

    F_likelihood((i - 1)*2+1 : (i-1)*2+2)      = F_likelihood_i;
    F_likelihood_mu_Y((i - 1)*2+1 : (i-1)*2+2) = F_likelihood_mu_Y_i;
end
%
drawBestTour(I_path , graph, J);

%+++++++++++++++++++++
Q_add = zeros(200,200);
Q_iminus{100} =pre_Q{100};
for i=1:N
    Q_add(((i-1)*2 + 1):((i-1)*2 + 2) , ((i-1)*2 + 1):((i-1)*2 + 2)) = Q_iminus{i};     % 
end


return
end



%%

function [V,Q] = path_prior_cov(N, F,pre_Q,R)
global Q_iminus
global PI
global Q_cov
global M_V


SIGMA_iminus         = zeros(2,2);   


T= diag([1.0 1.0].^2);
V = 0.0;
PI = 0.0;


for i=1:N-1

    % Jacobian, computed from the second node
    F_i     = F((i - 1)*2+1 : (i-1)*2+2);

    
    SIGMA_i_iminus = pre_Q{i};
    SIGMA_i_i      = SIGMA_i_iminus + SIGMA_iminus + T;

    SIGMA_i         = SIGMA_i_i - SIGMA_i_i*F_i'*invChol_mex((F_i*SIGMA_i_i*F_i' + R + eye(size(R))*1.0e-10))*F_i*SIGMA_i_i;
    SIGMA_i_i      = SIGMA_i_iminus + SIGMA_iminus + T ;


     
%   

    SIGMA_iminus = SIGMA_i;
    

    Q{i}        = SIGMA_i_i;        
   
    Q_iminus{i} = SIGMA_iminus;

    
    V = V + (F_i*SIGMA_i_i*F_i' + R);
    
    PI = PI + (F_i*SIGMA_i_iminus*F_i' + R);


    M_V(i) = (F_i*SIGMA_i_i*F_i' + R);

    
end
F_i     = F(99 : 100);
M_V(100) = (F_i*SIGMA_i_i*F_i' + R);

return
end



%%
function [X_path,I_path] = select(N,bar_X,Q,graph)
dimension = 2;

I_start = 1;
I_final = 1;
% 


I_node       = [[I_start:1:graph.n] ];
I_i          = I_start;                % the first node
I_node(I_i)  = [];               % list of the candidate node 

%
I_path(1) = I_start;
for i=1:N-1
    X_i         = [(graph.dist_x(I_i,I_node))' , (graph.dist_y(I_i,I_node))']; % distance to preceding node
    bar_X_i     = bar_X((i - 1)*dimension + 1 : (i - 1)*dimension + 2);        % expected distance, always greater than zero




    
    % select the node with the maximum likelihood
    N_i = length(I_node) ; 

    SIGMA_i_i   = Q{i};     % covariance of the distance
    L_i = Likelihood_of_node(N_i,X_i,bar_X_i,SIGMA_i_i);
    
    

    
    [L_i_min, I_min] = min(L_i);

    % update the cnadidate node
    I_i           = I_node(I_min);
    I_node(I_min) = [];

    % store the node in the path
    I_path(i)     = I_i;


end

data_L_I = [L_i',I_node',X_i];

I_path = [I_start,I_path,I_final];
X_path  = [];
for i=1:N
        X_path = [X_path , [graph.dist_x(I_path(i),I_path(i+1)) , graph.dist_y(I_path(i),I_path(i+1))]];
end


return
end






