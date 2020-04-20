%Constraint Function
   function [c, ceq] = Constraints_gradient_scaled(x)
       [cost,dT_chip,shear] = res_network_gradient_scaled(x(1),x(2),x(3),x(4));
   c = [(x(3) - x(2)) + 0.0002;
       (x(4) - x(3));
       shear - 37e6;
       dT_chip - 60];
   
   ceq = [];
   end