function [grad] = adjoint_gradient (adjoint_DSMC_conditions)

% This function computes the gradients of objective function with
% respect to a variable. The input variable adjoint_DSMC_conditions is a
% structure that includes 

% (a) gradient_type: this specifies the variable that the gradient is
% computed with respect to; 1 stands for gradient with respect to initial
% variance of particles velocities; 

% (b) bc: boundary condition. bc=1 stands for periodic
% boundary condition, bc=2 stands for reflecting boundary condition, bc=3
% stands for thermal boundary condition; 

% (c) VF: the velocities of particles at t=tmax; 

% (d) VI: the velocities of particles at t=0; 

% (e) XF: the locations of particles at t=tmax; 

% (f) bd: boundary domain. The boundary domain is specified as a 2 times
% DIMS matrix, where DIMS is the dimension of the spacial domain. The ith
% column represents the range of the spacial domain in the ith dimension;

% (g) tmax: the final time of the particles collisions;

% (h) collision_history: data of all particles collision; 

% (i) temperature: the initial variance of particles
% velocities. The output gradJ is the adjoint gradient.

% (j) inside_bd: indicator of whether each particle is inside boundary
% domain in each direction. This variable only applies to
% reflecting/thermal boundary conditions (i.e. bc=2 or bc=3)

% (k) V_ratio: ratio of velocities of post-updated velocities to
% pre-updated velocities; This variable only applies to
% thermal boundary conditions (i.e. bc=3)

% (l) pseudo_dt: the "delta t" for thermal boundary condition. This
% variable only applies to thermal boundary conditions （i.e. bc=3)

    dims = size(adjoint_DSMC_conditions.XF,2); %dimension of spacial domain
    r = adjoint_DSMC_conditions.velocity_function;
    drdv = adjoint_DSMC_conditions.velocity_function_derivative;
    f = adjoint_DSMC_conditions.measure_function;
    dfdx = adjoint_DSMC_conditions.measure_function_derivative;
    
    VF = adjoint_DSMC_conditions.VF;
    XF = adjoint_DSMC_conditions.XF;
    
    %set up initial adjoint variables
    switch dims
        case 1
            gamma = -drdv(VF).*f(XF);
            y = -r(VF).*dfdx(XF);
            
        case 2
            gamma = -drdv(VF).*f(XF(:,1),XF(:,2));
            y = -r(VF).*[dfdx(XF(:,1),XF(:,2)) dfdx(XF(:,2),XF(:,1))];
        
        case 3
            gamma = -drdv(VF).*f(XF(:,1),XF(:,2),XF(:,3));
            y = -r(VF).*[dfdx(XF(:,1),XF(:,2),XF(:,3)) dfdx(XF(:,2),XF(:,3),XF(:,1)) dfdx(XF(:,3),XF(:,1),XF(:,2))];
    
    end
    
    %compute gamma and y backwardly
    adjoint_DSMC_conditions.gamma = gamma;
    adjoint_DSMC_conditions.y = y;
    %{
    propagation_inputs.dt = adjoint_DSMC_conditions.dt;
    propagation_inputs.collision_history = adjoint_DSMC_conditions.collision_history;
    propagation_inputs.bc = adjoint_DSMC_conditions.bc;
    
    switch propagation_inputs.bc 
        case 2
            propagation_inputs.inside_bd = adjoint_DSMC_conditions.inside_bd;
        case 3
            propagation_inputs.inside_bd = adjoint_DSMC_conditions.inside_bd;
            propagation_inputs.V_ratio = adjoint_DSMC_conditions.V_ratio;
            propagation_inputs.pseudo_dt = adjoint_DSMC_conditions.pseudo_dt;
    end
    %}
    [gamma,y] = back_propagation(adjoint_DSMC_conditions);
    
    %compute gradients with respect to initial variance of particles velocities
    gradient_type = adjoint_DSMC_conditions.gradient_type;
    if (gradient_type==1) 
        grad = -sum(gamma.*(adjoint_DSMC_conditions.VI./(2*adjoint_DSMC_conditions.temperature)))/size(VF,1);
    end
    
    
end