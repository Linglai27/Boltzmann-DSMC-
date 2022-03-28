%% set up basic conditions for DSMC simulations and compute objective gradients
clear;

tic
nrun = 2; %number of runs of adjoint DSMC methods and centered difference methods
perturb = 0.05; %perturbation in centered difference methods
start_index = 0; %initial random seed number

adjoint_grad = zeros(nrun,3); %adjoint DSMC gradient record
finite_dif_grad = zeros(nrun,1); %centered difference gradient record

DSMC_inputs.N = 1e6; %number of particles involved in simulations
DSMC_inputs.nsteps = 1; %number of time steps
DSMC_inputs.temperature = [0.3 0.5 1]; %initial variance of particles' velocities
DSMC_inputs.dt = 0.1; %time step size

L = 2; %length of spacial domain in each direction
DSMC_inputs.bd = [0 0; L L]; %spacial/boundary domain
DSMC_inputs.uniform = 1; %indicate whether the initial location distribution of particles is uniform
DSMC_inputs.bc = 1; %the type of boundary condition(1=periodic, 2=reflecting, 3=thermal)

%set up boundary temperature if thermal condition is used
if (DSMC_inputs.bc==3)
    if (size(DSMC_inputs.bd,2)==1) %1D boundary domain
        %left boundary temperature for 1D thermal case
        DSMC_inputs.left_boundary_temperature = [0.3 0.4 0.5];

        %right boundary temperature for 1D thermal case
        DSMC_inputs.right_boundary_temperature = [0.5 0.7 0.9];
    else %2D/3D boundary domain
        DSMC_inputs.boundary_temperature = [0.5 0.7 0.9];
    end
end

DSMC_inputs.alpha = 0; %type of particle collision, alpha = 0 means Maxwellian
DSMC_inputs.mu = 1; DSMC_inputs.eps = 1; DSMC_inputs.rho = 1;
DSMC_inputs.Nc = 10; %number of cells in each spacial domain direction

DSMC_inputs.velocity_function = @(v) sum(v.^4,2); %velocity function r
DSMC_inputs.velocity_function_derivative = @(v) 4*v.^3; %velocity function derivative

%measure domain function and its derivative
[DSMC_inputs.measure_function, DSMC_inputs.measure_function_derivative]...
    = measure_function(size(DSMC_inputs.bd,2), [0 0], 0.01, [L/DSMC_inputs.Nc L/DSMC_inputs.Nc]*3); 
DSMC_inputs.gradient_type = 1; %type of gradient, 1 stands for gradient with respect to initial temperature of particles

%% compute adjoint gradients
parfor i = 1:nrun
    adjoint_grad(i,:) = adjoint_gradient(NH_solver(DSMC_inputs,start_index+i));
end

DSMC_inputs.temperature = DSMC_inputs.temperature + [0 0 perturb];

%% compute finite difference gradient
parfor i = 1:nrun
    finite_dif_grad(i,:) = objective_func(NH_solver(DSMC_inputs,start_index+i));
end

DSMC_inputs.temperature = DSMC_inputs.temperature - 2*[0 0 perturb];

parfor i = 1:nrun
    finite_dif_grad(i,:) = (finite_dif_grad(i,:) - objective_func(NH_solver(DSMC_inputs,start_index+i)))/(2*perturb);
end

%% print means and variance of gradients

fprintf('The finite difference gradient mean is %f \n', vpa(mean(finite_dif_grad),6));
m = mean(adjoint_grad,1);
fprintf('The adjoint gradient mean is [%f %f %f] \n', vpa(m(1),6), vpa(m(2),6), vpa(m(3),6));

fprintf('The finite difference gradient variance is %f \n', vpa(var(finite_dif_grad),10));
v = var(adjoint_grad,1);
fprintf('The adjoint gradient variance is [%f %f %f] \n', vpa(v(1),10), vpa(v(2),10), vpa(v(3),10));
toc