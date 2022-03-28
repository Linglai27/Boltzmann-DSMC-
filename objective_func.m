function J = objective_func(DSMC_outputs)

    % This function computes the values of objective function which will be
    % used to compute the finite difference gradients. The input
    % DSMC_outputs is a structure including 
    
    % (a) VF: the final velocities of the particles; 
    
    % (b) XF: the final locations of the particles; 
    
    % (c) r: the velocity function;
    
    % (d) I: the measure domain indicator function.
    
    % The output J is the value of the objective function.
    
    VF = DSMC_outputs.VF;
    XF = DSMC_outputs.XF;
    r = DSMC_outputs.velocity_function;
    I = DSMC_outputs.measure_function;
    
    dims = size(XF,2);
    
    switch dims
        case 1
        J = mean(r(VF).*I(XF));
        case 2
        J = mean(r(VF).*I(XF(:,1),XF(:,2)));
        case 3
        J = mean(r(VF).*I(XF(:,1),XF(:,2),XF(:,3)));
    end
end