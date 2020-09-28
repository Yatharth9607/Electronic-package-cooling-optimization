%Application of GA
%clear
%clc
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

% plot of the effect of design variables on fitness function
% Number of solder balls
var = linspace(x(1) * 0.8, x(1) * 1.2, 10);
for i = 1:10
    [cost_nb(i), dT_chip_nb(i)] = res_network_gradient(var(i), x(2), x(3), x(4));    
end

figure(1)
yyaxis left
plot(var, cost_nb)
ylabel('Cost of the module')

hold on
yyaxis right
plot(var, dT_chip_nb)
ylabel(['Temperature rise in chip (' char(176) 'C)'])

xlabel('Number of solder balls')

% Pitch of solder array
var = linspace(x(2) * 0.8, x(2) * 1.2, 10);
for i = 1:10
    [cost_nb(i), dT_chip_nb(i)] = res_network_gradient(x(1), var(i), x(3), x(4));    
end

figure(2)
yyaxis left
plot(var, cost_nb)
ylabel('Cost of the module')

hold on
yyaxis right
plot(var, dT_chip_nb)
ylabel(['Temperature rise in chip (' char(176) 'C)'])

xlabel('Pitch of the solder ball array (m)')

% Diameter of solder balls
var = linspace(x(3) * 0.8, x(3) * 1.2, 10);
for i = 1:10
    [cost_nb(i), dT_chip_nb(i)] = res_network_gradient(x(1), x(2), var(i), x(4));    
end

figure(3)
yyaxis left
plot(var, cost_nb)
ylabel('Cost of the module')

hold on
yyaxis right
plot(var, dT_chip_nb)
ylabel(['Temperature rise in chip (' char(176) 'C)'])

xlabel('Diameter of solder balls (m)')

% Height of solder balls
var = linspace(x(4) * 0.8, x(4) * 1.2, 10);
for i = 1:10
    [cost_nb(i), dT_chip_nb(i)] = res_network_gradient(x(1), x(2), x(3), var(i));    
end

figure(4)
yyaxis left
plot(var, cost_nb)
ylabel('Cost of the module')

hold on
yyaxis right
plot(var, dT_chip_nb)
ylabel(['Temperature rise in chip (' char(176) 'C)'])

xlabel('Height solder balls (m)')