function y = fmincon_scaled_objective(x)
%     Nb = x(1)*1500;
%     p = x(2)*0.0002769;
%     db = x(3)*0.0005073;
%     hb = x(4)*5.75044e-05;
    Nb = x(1);
    p = x(2);
    db = x(3);
    hb = x(4);
    [cost,dT_chip] = res_network_gradient_scaled(Nb,p,db,hb);

%    y = dT_chip;      % Single Objective - Temperature change of Chip
    a = 1/3;
    b = 2/3;
    y = a*cost + b*dT_chip;     % Multi-objective
end