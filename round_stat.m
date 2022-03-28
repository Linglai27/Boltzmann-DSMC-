function out = round_stat(x)
%ROUND_STAT  Statistical rounding
% out = [x]+1   with probability x-[x],  [x]=floor(x)
% out = [x]     with probability 1-(x-[x])
% so that out is a vector of integers and mean(out)=mean(x)
r=rand(size(x));
fx=floor(x);
out = (fx+1).*(r<=x-fx)  + fx.*(r>x-fx);
end

