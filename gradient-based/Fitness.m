%Fitness Function
   function y = Fitness(x)
    %y = 100 * (x(1)^2 - x(2)) ^2 + (1 - x(1))^2;
    Nb = x(1);
    p = x(2);
    db = x(3);
    [total_cost,dT] = res_network(Nb,p,db);
    a = 1/3;
    b = 2/3;
    y = a*total_cost + b*dT;      %Maximization for a>0 and b<0
   end