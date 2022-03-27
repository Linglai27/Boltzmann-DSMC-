function DSMC = NH_solver(DSMC, index)

% This function simulate gas particles under spacially inhomogeneous
% Boltzmann equation with different boundary conditions (periodic,
% reflecting and thermal) in 1,2 or 3 dimensions rectangular spacial
% domain. The DSMC is a structure variable which includes the 
% following name needed are 

% (a) bc: boundary condition. bc=1 stands for periodic
% boundary condition, bc=2 stands for reflecting boundary condition, bc=3
% stands for thermal boundary condition. 

% (b) bd: boundary domain. The boundary domain is specified as a 2 times
% DIMS matrix, where DIMS is the dimension of the spacial domain. The ith
% column represents the range of the spacial domain in the ith dimension;

% (c) N: the number of particles involved in the 
% simulations; 

% (d) nsteps: number of discretized time steps; 

% (e) dt: the size of the discretized time step; 

% (f) temperature: the initial temperature (variance of particles) in each 
% dimension; 

% (g) alpha: the type of collision rule for the gas. Alpha = 0 stands for
% Maxwellian gas and Alpha=1 stands for VHS gas; 

% (h) mu: a parameter impacting proportion of colliding particles when
% alpha = 0;

% (i) eps: a parameter impacting proportion of colliding particles when
% alpha = 0;

% (j) rho: a parameter determining proportion of colliding particles when
% alpha = 1;

% (h) index: the random seed used; 

% (i) Nc: the number of cells in each direction; 

% (j) uniform: determine whether the initial location distribution of
% particles is uniform. uniform=1 stands for uniform distribution,
% uniform=0 stands for normal distribution;

% (k) left_boundary_temperature: the temperature for particles out of left
% boundary in 1D thermal case;

% (l) right_boundary_temperature: the temperature for particles out of
% right boundary in 2D/3D thermal case;

% (m) boundary_temperature: the temperature for particles out of boundary
% in 2D/3D thermal cases.

% (n) VF: velocities of particles at t=tmax=nsteps*dt; 

% (o) XF: locations of particles at t=tmax; 

% (p) collision_history: data of all particles collision; 

% (q) VI: velocities of particles at t=0; 

% (r) XI: locations of particles at t=0;

% (s) inside_bd: indicator of whether each particle is inside boundary
% domain in each direction. This variable only applies to
% reflecting/thermal boundary conditions (i.e. bc=2 or bc=3)

% (t) V_ratio: ratio of velocities of post-updated velocities to
% pre-updated velocities; This variable only applies to
% thermal boundary conditions (i.e. bc=3)

% (u) pseudo_dt: the "delta t" for thermal boundary condition. This
% variable only applies to thermal boundary conditions （i.e. bc=3)

% The second input index specifies the random seed number.
    
    rng default;
    rng(index);
    
    T0 = DSMC.temperature;
    bd = DSMC.bd;
    N = DSMC.N;
    
    V = standard_randn(N,3).*sqrt(T0);
    DSMC.VI = V;
    uniform = DSMC.uniform; 
    dims = size(bd,2);

    %set up initial location distribution
    if (uniform)
        X = bd(1,:) + 1/3*(bd(2,:)-bd(1,:)).*rand(N,dims);
    else %otherwise normal initial location distribution
        X = (bd(1,:)+bd(2,:))/2+(bd(2,:)-bd(1,:))/5*standard_randn(N,dims);
    end
    
    DSMC.XI = X;
    
    Nc = DSMC.Nc; %number of cell in each dimension
    dt = DSMC.dt;
    nsteps = DSMC.nsteps;
    alpha = DSMC.alpha;
    
    mu = DSMC.mu; eps = DSMC.eps; rho = DSMC.rho;

    collision_history = cell(nsteps,3);
    particle_index = 1:N;
    
    bc = DSMC.bc;
    
    switch bc
        case 2 
            DSMC.inside_bd = zeros(N,dims,nsteps);
        case 3
            DSMC.inside_bd = zeros(N,dims,nsteps);
            DSMC.pseudo_dt = zeros(N,dims,nsteps);
            DSMC.V_ratio = zeros(N,dims,nsteps);
    end
    
    for i = 1:DSMC.nsteps
        for j = 1:Nc^dims
            [ind] = compute_index(j,Nc,dims);
            selected = ones(N,1);
            for k = 1:dims
                L = bd(2,k)-bd(1,k);
                selected = selected&(X(:,k)/L>((ind(k)-1)/Nc))&(X(:,k)/L<=(ind(k)/Nc));
            end

            if (sum(selected)>1)
              [~,V(selected,:),~,u_dir,sigma,collide_index]= H_solver(V(selected,:),dt,mu,rho,eps,alpha,particle_index(selected));
            end
            
            collision_history{i,1} = [collision_history{i,1}; u_dir];
            collision_history{i,2} = [collision_history{i,2}; sigma];
            
            tmp = collision_history{i,3};
            collision_history{i,3} = [tmp(1:length(tmp)/2) collide_index(1:length(collide_index)/2)...
                tmp(length(tmp)/2+1:length(tmp)) collide_index(length(collide_index)/2+1: length(collide_index))];
            
        end
        
        X = X + V(:,1:dims).*dt;
        
        bc_inputs.N = N;
        bc_inputs.V = V;
        bc_inputs.X = X;
        bc_inputs.bd = bd;
        bc_inputs.bc = bc;
        
        if (bc==3)
            if (dims==1)
                bc_inputs.left_boundary_temperature = DSMC.left_boundary_temperature;
                bc_inputs.right_boundary_temperature = DSMC.right_boundary_temperature;
            else
                bc_inputs.boundary_temperature =  DSMC.boundary_temperature;
            end
        end
        
        bc_outputs = boundary_condition(bc_inputs);
        
        V = bc_outputs.V;
        X = bc_outputs.X;
        
        switch bc
            case 2
                DSMC.inside_bd(1:N,1:dims,i) = bc_outputs.inside_bd;
            case 3
                DSMC.inside_bd(1:N,1:dims,i) = bc_outputs.inside_bd;
                DSMC.pseudo_dt(1:N,1:dims,i) = bc_outputs.pseudo_dt + dt;
                DSMC.V_ratio(1:N,1:dims,i) = bc_outputs.V_ratio;
        end
    end

    DSMC.collision_history = collision_history;
    DSMC.VF = V;
    DSMC.XF = X;
end

