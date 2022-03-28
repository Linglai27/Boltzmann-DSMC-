function [gamma, y] = back_propagation(propagation_inputs)

% This function computes the gamma and y backwardly for different boundary
% conditions. The input propagation_inputs is a structure including

% (a) bc: boundary condition. 1 stands for periodic boundary condition, 2
% stands for reflecting boundary condition, 3 stands for thermal boundary
% condition. (b) gamma: gammas for particles at t=tmax; 

% (b) y: ys for particles at t=tmax; 

% (c) tmax: the final time of particles motions; 

% (d) nsteps: number of discretized time steps; 

% (e) collision_history: data of all particles collision; 

% (f) inside_bd: indicator of whether each particle is inside boundary
% domain in each direction. This variable only applies to
% reflecting/thermal boundary conditions (i.e. bc=2 or bc=3)

% (g) V_ratio: ratio of velocities of post-updated velocities to
% pre-updated velocities; This variable only applies to
% thermal boundary conditions (i.e. bc=3)

% (h) pseudo_dt: the "delta t" for thermal boundary condition. This
% variable only applies to thermal boundary conditions （i.e. bc=3)

% The outputs are 

% (a) gamma: gammas for particles at t=0; 

% (b) y: ys for particles at t=0.
    
    gamma = propagation_inputs.gamma;
    y = propagation_inputs.y;
    N = size(gamma,1); %number of particles involved in the simulation
    dims = size(y,2); %dimension of spacial domain
    dt = propagation_inputs.dt; 
    bc = propagation_inputs.bc;
    nsteps = size(propagation_inputs.collision_history,1);
    
    %start backward propagation computation
    for i = nsteps:-1:1
        uncollided = ones(N,1);
        
        collided = propagation_inputs.collision_history{i,3};
            
        uncollided(collided') = 0;
        
        switch bc
            case 1 %process gamma for colliding particles under periodic boundary condition
            
                gamma(collided,:) = COLLIDEgamma(gamma(collided,:)+[dt.*ones(size(collided,2),dims).*y(collided,:) zeros(size(collided,2),3-dims)],...
                    propagation_inputs.collision_history{i,2},propagation_inputs.collision_history{i,1});
                gamma(uncollided==1,:) = gamma(uncollided==1,:)+[dt.*ones(sum(uncollided),dims).*y(uncollided==1,:) zeros(sum(uncollided),3-dims)];

            case 2 %process gamma for colliding particles under reflecting boundary condition
            
                reflect_mat = 2*propagation_inputs.inside_bd(collided,:,i)-1;
                gamma(collided,:) = COLLIDEgamma([reflect_mat ones(size(collided,2),3-dims)].*(gamma(collided,:)+...
                    [dt.*ones(size(collided,2),dims).*y(collided,:) zeros(size(collided,2),3-dims)]),...
                    propagation_inputs.collision_history{i,2},propagation_inputs.collision_history{i,1});

                reflect_mat = 2*propagation_inputs.inside_bd(uncollided==1,:,i)-1;
                gamma(uncollided==1,:) = [reflect_mat ones(sum(uncollided),3-dims)].*(gamma(uncollided==1,:)+...
                    [dt.*ones(sum(uncollided),dims).*y(uncollided==1,:) zeros(sum(uncollided),3-dims)]);

                y = (propagation_inputs.inside_bd(:,1:dims,i)).*y;

            case 3 %process gamma for colliding particles under thermal boundary condition
                
                inside_bd = propagation_inputs.inside_bd(:,:,i);
                pseudo_dt = propagation_inputs.pseudo_dt(:,:,i);
                V_ratio = propagation_inputs.V_ratio(:,:,i);
                
                gamma(collided,:) = COLLIDEgamma([inside_bd(collided,1:dims)==0 ones(size(collided,2),3-dims)].*gamma(collided,:)+...
                        [pseudo_dt(collided,:).*V_ratio(collided,:).*y(collided,1:dims)  zeros(size(collided,2),3-dims)],...
                        propagation_inputs.collision_history{i,2},propagation_inputs.collision_history{i,1});

                gamma(uncollided==1,:) = [inside_bd(uncollided==1,1:dims)==0 ones(sum(uncollided),3-dims)].*gamma(uncollided==1,:)+...
                   [pseudo_dt(uncollided==1,:).*V_ratio(uncollided==1,:).*y(uncollided==1,1:dims)  zeros(sum(uncollided),3-dims)];

                y = V_ratio(:,:).*y;
        end
        
    end
end
