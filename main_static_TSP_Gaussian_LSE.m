function main_static_TSP_Gaussian_LSE

% The matlab code for mixed integere TSP construction.  The continuous-based
% expression and estmate-correct mechnism is used. 
%

% The author: 
% Liqiang Hou (Dr.Who)
% Research Prof., Ph.D
% Shanghai Jioatong University

% All rights reserved

clear all;
close all;
global M_T


addpath(genpath('branch-03-15'));
graph = construct_graph_kroA100();


% % 21282, the actual optimal
% I_Path_21282 = [1	47	93	28	67	58	61	51	87	25	81	69	64	40	54	2	44	50	73	68	85	82	95	13	76	33	37	5	52	78	96	39	30	48	100	41	71	14	3	43	46	29	34	83	55	7	9	57	20	12	27	86	35	62	60	77	23	98	91	45	32	11	15	17	59	74	21	72	10	84	36	99	38	24	18	79	53	88	16	94	22	70	66	26	65	4	97	56	80	31	89	42	8	92	75	19	90	49	6	63	1];
% drawBestTour(I_Path_21282, graph, 21285.44);


tic

[bar_X,W,R,M_T] = Init_X_Set(graph);

J = main_fun(bar_X ,W,R, graph);

toc
return
end


%%

function J = main_fun(bar_X, W ,R, graph)
global M_T
global Q_iminus

global Q_add
global bar_J
global R
global F_likelihood


delta_J = 1.0e3;
N       = graph.n; 
M_T = M_T*0.5;

%

for i=1:N
    Q_iminus{i} = zeros(2,2);     % 
end




% Initilialize Jacobian

X_path = bar_X;
P_X    = W;

[F]                 = path_Jacobian(N, bar_X);
[J,X_path,V,I_path] = cost_fun_tsp_estimation(bar_X, F,P_X,R, graph);


pre_bar_J = bar_J;
pre_J     = J;
pre_V     = V; 

%
bar_X = bar_X + ones(size(bar_X))*30.0;     % 



while delta_J > 1.0e-3

    % update the Jacobian
    [F_apth]                    = path_Jacobian(N, X_path);
    [J,X_path,V,I_path]    = cost_fun_tsp_estimation(bar_X,F_apth,P_X,R, graph);

    H    = F_likelihood';

    P_XY = P_X*H;
    P_YX = (P_XY)';
    P_Y  = H'*(P_X + M_T)*H + R;
  


    % conditional predict
    P_X   = P_X - P_XY/P_Y*P_YX;

    bar_X = bar_X + P_XY/P_Y*( (0.5*(bar_J - J)/V*(bar_J - J) + J) ...
          - (0.5*(pre_bar_J - pre_J)/V*(pre_bar_J - pre_J) + pre_J));     

    
    % fusion
    % 
     P_X     = pinv( ...
               pinv((M_T + eye(size(W))*1.0e-8),   1.0e-10) ...
             + pinv((Q_add + eye(size(W))*1.0e-8), 1.0e-10) ...
             + pinv((P_X  + eye(size(W))*1.0e-8) , 1.0e-10) ...
             + eye(size(W))*1.0e-8             ,   1.0e-10);


     % 

    delta_J = abs(J - pre_J);
    pre_bar_J = bar_J;

    pre_J = J; 
    pre_V = V;
end


% draw best tour

drawBestTour(I_path , graph, J);

return
end




