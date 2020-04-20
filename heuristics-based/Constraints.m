%Constraint Function
   function [c, ceq] = Constraints(x)
   c = (x(3) - x(2)) + 2;
   ceq = [];
   end