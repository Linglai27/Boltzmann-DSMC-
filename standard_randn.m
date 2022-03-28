function real_normal = standard_randn(Nx,Ny)
    
    % This function adjusts random numbers from randn() so that each column
    % has mean 0 and variance 1 within machine precision

	real_normal = randn(Nx,Ny);
	real_normal = real_normal-mean(real_normal);
	real_normal = real_normal./std(real_normal);
end