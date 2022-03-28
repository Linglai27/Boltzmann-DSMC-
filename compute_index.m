function ind = compute_index(N,Nc,dims)
    if (dims==1)
        ind = N;
    else
        tmp = ceil(N/Nc^(dims-1));

        ind = [tmp; compute_index(N - (tmp-1)*Nc^(dims-1),Nc,dims-1)];
    end
end