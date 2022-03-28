function bc_outputs  = boundary_condition(bc_inputs)

    % This function updates velocities and locations of particles based on
    % bd conditions and bd domains. The input bc_inputs is a input
    % structure which includes 
    
    % (a) bc: the type of bd condition used. 1 stands for periodic bd
    % condition; 2 stands for reflecting bd condition; 3 stands for
    % thermal bd condition. For thermal case, bc{2,1} includes the boundary
    % temperature;
    
    % (b) X: the current locations of particles; 
    
    % (c) V: the current velocities of particles; 
    
    % (d) bd: boundary domain. This parameter specifies the bd domain. Each
    % row represents the boundary domain on the specified dimension. For
    % instance, [1,2;3,4] represents the domain [1,3] \times [2,4]. 
    
    % (e) left_boundary_temperature: the temperature for particles out of
    % left boundary in 1D thermal case;

    % (f) right_boundary_temperature: the temperature for particles out of
    % right boundary in 1D thermal case;

    % (g) boundary_temperature: the temperature for particles out of
    % boundary in 2D/3D thermal cases.
    
    % The outputs are 
    
    % (a) X: the updated locations of particles based on boundary
    % condition; 
    
    % (b) V: the updated velocities of particles based on
    % boundary condition; 
    
    % (c) inside_bd: indicator of whether each particle is inside boundary
    % domain in each direction. This variable only applies to
    % reflecting/thermal boundary conditions (i.e. bc=2 or bc=3)

    % (d) V_ratio: ratio of velocities of post-updated velocities to
    % pre-updated velocities; This variable only applies to
    % thermal boundary conditions (i.e. bc=3)

    % (e) pseudo_dt: the "delta t" for thermal boundary condition. This
    % variable only applies to thermal boundary conditions （i.e. bc=3)

    bd = bc_inputs.bd;
    dims = size(bd,2); %dimension of spacial domain
    bc = bc_inputs.bc;
    X = bc_inputs.X;
    V = bc_inputs.V;
    N = bc_inputs.N;
    
    switch bc
        
        case 1 %periodic bd condition
        
        for i = 1:dims
            X(:,i) = bd(1,i)+mod(X(:,i),bd(2,i)-bd(1,i));
        end
            
        case 2
        inside_bd = zeros(size(V,1),dims);
        
        for i = 1:dims
            out_of_bd_1 = (X(:,i)<bd(1,i)); 
            out_of_bd_2 = (X(:,i)>bd(2,i)); 

            V((out_of_bd_1|out_of_bd_2),i) = -V((out_of_bd_1|out_of_bd_2),i);
            X(out_of_bd_1,i) = 2*bd(1,i) - X(out_of_bd_1,i);
            X(out_of_bd_2,i) = 2*bd(2,i) - X(out_of_bd_2,i);
            
            inside_bd(:,i) = ~(out_of_bd_1|out_of_bd_2); %1 if inside bd, 0 otherwise
        end
        
        bc_outputs.inside_bd = inside_bd;
    
        case 3
        
        if (dims==1)
            L = bd(1,1); R = bd(2,1);
            out_of_bd_1 = (X<L);
            out_of_bd_2 = (X>R);
            
            pseudo_dt = zeros(N,1); %the new delta t for particles out of bd
            pseudo_dt(out_of_bd_1,1) = (L-X(out_of_bd_1))./V(out_of_bd_1,1);
            pseudo_dt(out_of_bd_2,1) = (R-X(out_of_bd_2))./V(out_of_bd_2,1);
            
            g = standard_randn(sum(out_of_bd_1,1),3).*sqrt(bc_inputs.left_boundary_temperature);
            h = standard_randn(sum(out_of_bd_2,1),3).*sqrt(bc_inputs.right_boundary_temperature);
            g(:,1) = abs(g(:,1));
            h(:,1) = -abs(h(:,1));
            
            V_ratio = ones(N,1);
            V_ratio(out_of_bd_1,1) = g(:,1)./V(out_of_bd_1,1);
            V_ratio(out_of_bd_2,1) = h(:,1)./V(out_of_bd_2,1);
            
            V(out_of_bd_1,:) = g;
            V(out_of_bd_2,:) = h;
            V(out_of_bd_1,1) = abs(V(out_of_bd_1,1));
            V(out_of_bd_2,1) = -abs(V(out_of_bd_2,1));
            
            X(out_of_bd_1,1) = L + pseudo_dt(out_of_bd_1,1).*V(out_of_bd_1,1);
            X(out_of_bd_2,1) = R + pseudo_dt(out_of_bd_2,1).*V(out_of_bd_2,1);
            
            %2 means out of right bd, 1 means out of left bd, 0 means inside the bd
            %domain
            inside_bd = ones(size(V,1),1);
            inside_bd(:,1) = 2*out_of_bd_2+out_of_bd_1; 

        else
            
            inside_bd = zeros(N,dims);
            pseudo_dt = zeros(N,dims); %the new delta t for particles out of bd
            g = standard_randn(N,3).*sqrt(bc_inputs.boundary_temperature);
            V_ratio = ones(N,dims);
            
            for i = 1:dims
                out_of_bd_1 = X(:,i)<bd(1,i); %identify particles out of left/lower/bottom boundary
                out_of_bd_2 = X(:,i)>bd(2,i); %identify particles out of right/upper/top boundary
                
                %0 if inside bd, 1 if out of left/lower/bottom
                %bd, 2 if out of right/upper/top bd
                inside_bd(:,i) = 2*out_of_bd_2+out_of_bd_1; 
                pseudo_dt(out_of_bd_1,i) = (bd(1,i)-X(out_of_bd_1,i))./V(out_of_bd_1,i);
                pseudo_dt(out_of_bd_2,i) = (bd(2,i)-X(out_of_bd_2,i))./V(out_of_bd_2,i);
                V_ratio(out_of_bd_1,i) = abs(g(out_of_bd_1,i))./V(out_of_bd_1,i);
                V_ratio(out_of_bd_2,i) = -abs(g(out_of_bd_2,i))./V(out_of_bd_2,i);
                V(out_of_bd_1,i) = abs(g(out_of_bd_1,i));
                V(out_of_bd_2,i) = -abs(g(out_of_bd_2,i));
                X(out_of_bd_1,i) = bd(1,i) + V(out_of_bd_1,i).*pseudo_dt(out_of_bd_1,i);
                X(out_of_bd_2,i) = bd(2,i) + V(out_of_bd_2,i).*pseudo_dt(out_of_bd_2,i);
            end
            
        end
        
        bc_outputs.inside_bd = inside_bd;
        bc_outputs.pseudo_dt = pseudo_dt;
        bc_outputs.V_ratio = V_ratio;
    end
    
    bc_outputs.V = V;
    bc_outputs.X = X;
end