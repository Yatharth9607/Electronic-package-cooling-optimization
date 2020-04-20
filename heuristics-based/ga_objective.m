%Fitness Function
   function y = ga_objective(x)

    [total_cost,dT] = res_network_heuristic(x(1),x(2),x(3),x(4));
    a = 1/3;
    b = 2/3;
    y = a*total_cost + b*dT;      %Maximization for a>0 and b<0
   end