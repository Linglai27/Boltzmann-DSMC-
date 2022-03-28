function velocity_function(V)
    global r
    global drdv
    
    r = sum(V.^4,2);
    drdv = sum(4*V.^3,2);
end