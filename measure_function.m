function [f, dfdx] = measure_function(dims, center, std, width)

    %This function returns the measure function f. The inputs are (a) dims:
    %the dimension of the boundary domain; (b) center: the center of the
    %measure domain; (c) std: the standard deviation of the normal
    %probabilistic distribution function used ; (d) width: the width of the
    %measure window. The outputs are (a) f: the measure function; (b) the
    %(partial) derivative of f with respect to the first variable of f.

    
    switch dims
        case 1
        f = @(x) -normcdf(x,center+width/2,std) + normcdf(x,center-width/2,std);
        dfdx = @(x) -normpdf(x,center+width/2,std) + normpdf(x,center-width/2,std);
    
        case 2
        f = @(x1,x2) -normcdf(x1,center(1)+width(1)/2,std).*normcdf(x2,center(2)+width(2)/2,std) + normcdf(x1,center(1)-width(1)/2,std).*normcdf(x2,center(2)-width(2)/2,std);
        dfdx = @(x1,x2) -normpdf(x1,center(1)+width(1)/2,std).*normcdf(x2,center(2)+width(2)/2,std) + normpdf(x1,center(1)-width(1)/2,std).*normcdf(x2,center(2)-width(2)/2,std);
        
        case 3
        f = @(x1,x2,x3) -normcdf(x1,center(1)+width(1)/2,std).*normcdf(x2,center(2)+width(2)/2,std).*normcdf(x3,center(3)+width(3)/2,std) + normcdf(x1,center(1)-width(1)/2,std).*normcdf(x2,center(2)-width(2)/2,std).*normcdf(x3,center(3)-width(3)/2,std);
        dfdx = @(x1,x2,x3) -normpdf(x1,center(1)+width(1)/2,std).*normcdf(x2,center(2)+width(2)/2,std).*normcdf(x3,center(3)+width(3)/2,std) + normpdf(x1,center(1)-width(1)/2,std).*normcdf(x2,center(2)-width(2)/2,std).*normcdf(x3,center(3)-width(3)/2,std);
    end
end
