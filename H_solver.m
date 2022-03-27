function [obj,VF,VI,u_dir,sigma,collide_index] = H_solver(v_particle,dt,mu,rho,eps,alpha,particle_index)
    VI = v_particle;
    num_particle = size(v_particle,1);
        
    if (alpha ==1) %VHS
    p = 1;
    sigma = (2*max(sqrt(v_particle(:,1).^2+v_particle(:,2).^2+v_particle(:,3).^2)))^p;
    num_pair = round_stat(num_particle*sigma*rho*dt/2/eps);
    perm = randperm(num_particle,2*num_pair);
    perm1 = perm(1:num_pair);
    perm2 = perm(num_pair+1:2*num_pair);
    V1 = v_particle(perm(1:num_pair),:);
    V2 = v_particle(perm(num_pair+1:2*num_pair),:);
    v_select = V1 - V2;
    sigmaij_square = v_select(:,1).^2+v_select(:,2).^2+v_select(:,3).^2;
    sigmaij = sqrt(sigmaij_square).^p;
    collision_index = (sigma.*rand(num_pair,1)<sigmaij)';
    perm = [perm1(collision_index) perm2(collision_index)];
    else %Maxwellian
    num_pair = round_stat(num_particle*mu/eps*dt/2);
    perm = randperm(num_particle,2*num_pair);
    end
    
    [v_particle(perm,:), u_dir, sigma] = Collision(v_particle(perm,:));
    collide_index = particle_index(perm);
    
    VF = v_particle;
    obj = [sum(v_particle.^2)/num_particle  sum(v_particle.^2.^2)/num_particle]';
end