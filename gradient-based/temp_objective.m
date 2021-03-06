function y = temp_objective(x)
    Nb = x(1);
    p = x(2);
    db = x(3);
    hb = x(4);
    [cost,dT_chip] = res_network_gradient(Nb,p,db,hb);

%    y = dT_chip;      % Single Objective - Temperature change of Chip
    a = 0;
    b = 1;
    y = a*cost + b*dT_chip;     % Multi-objective
end