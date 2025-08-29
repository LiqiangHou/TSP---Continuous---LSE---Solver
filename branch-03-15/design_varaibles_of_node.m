function [mu_x,SIGMA_x,kappa] = design_varaibles_of_node(x,i_node,SIGMA_x_pre,graph,W)
    
    n_var        = 4;
    
    m_target     = x((i_node - 1)*n_var + 1  : (i_node - 1)*n_var + 2);         % dx,dy
    d_target     = x((i_node - 1)*n_var + 3  : (i_node - 1)*n_var + 4);         % std_x,std_y
    


 

    min_d_pos     = [graph.min_dx        graph.min_dy];
    max_d_pos     = [graph.max_dx        graph.max_dy];
    max_std_dr    = [graph.max_abs_dst_x  ,...
                     graph.max_abs_dst_y  ];      

    %
    mu_x     = min_d_pos + (max_d_pos - min_d_pos).*m_target;   % from [0  1] to [min_dr max_dr]
    
    SIGMA_x  = diag(d_target .* max_std_dr ) + W + SIGMA_x_pre;


return
end