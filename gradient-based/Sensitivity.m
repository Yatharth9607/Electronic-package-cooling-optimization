%Application of GA
clear
clc
clear; clc;
x0 = [511 0.0010 0.00087 0.00055]; % initial point [Nb,p,db,hb]
A = [];
b = [];
Aeq = [];
beq = [];
lb = [25 0.0008 0.0006 0.0002];   % Lower bounds
ub = [2500 0.0020 0.0020 0.0020];  % Upper bounds

ConstraintFunction = @Constraints_gradient;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

objective_function = @fmincon_objective; % objective function

[x,Fval] = fmincon(objective_function,x0,A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
[cost,dT_chip] = res_network_gradient(x(1),x(2),x(3),x(4));
% Jacobian of multi-objective function
x0 = x;
dx = 0.00001;
for i = 1:4
    
    a = x0(i) - dx;
    b = x0(i) + dx;
    
    xa = x0;
    xb = x0;
    
    xa(i) = a;
    xb(i) = b;
    
    [cost_a,dT_chip_a] = res_network_gradient(xa(1),xa(2),xa(3),xa(4));
    [cost_b,dT_chip_b] = res_network_gradient(xb(1),xb(2),xb(3),xb(4));
    J_cost(i) = (cost_b - cost_a)/(2*dx);
    J_dT_chip(i) = (dT_chip_b - dT_chip_a)/(2*dx);
    
    mult_obj_a = (1/3)*cost_a + (2/3)*dT_chip_a;
    mult_obj_b = (1/3)*cost_b + (2/3)*dT_chip_b;
    J_obj(i) = (mult_obj_b - mult_obj_a)/(2*dx);
    
end


J_norm_cost(1,:) = J_cost(:).*x0(:)/cost;
J_norm_dT_chip(1,:) = J_dT_chip(:).*x0(:)/dT_chip;
J_norm_obj(1,:) = J_obj(:).*x0(:)/Fval;

disp(J_norm_cost)
disp(J_norm_dT_chip)
disp(J_norm_obj)