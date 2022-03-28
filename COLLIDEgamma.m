% gamma-collisions
% This is effectively applying d V'_{k,i}/ dV_{k,i}, etc.
function gamma_out=COLLIDEgamma(gamma, omegas, u_dir)
    N=size(gamma,1);
    gamma_avg  = (gamma(1:N/2,:) + gamma(N/2+1:N,:))/2;
    gamma_diff = (gamma(1:N/2,:) - gamma(N/2+1:N,:))/2;
    
    Dg=sum(gamma_diff.*omegas,2).*u_dir;
    
    gamma_out=[gamma_avg; gamma_avg] + [Dg;-Dg];
end

