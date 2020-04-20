% % Pareto front - Normal Boundary Intersection
clear; clc;
N = 51;
x0 = [0,0]; % initial point [Nb,p,db]
A = [];
b = [];
Aeq = [];
beq = [];
lb = [-3 -3];   % Lower bounds
ub = [3 3];  % Upper bounds

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

objective_function_1 = @objective_1; % first objective function
objective_function_2 = @objective_2; % second objective function

% Shadow minima
[x_1,Fval_1] = fmincon(objective_function_1,x0,A,b,Aeq,beq,lb,ub);   % calling the optimizer
[J(1,1),J(1,2)] = example(x_1);
[x_2,Fval_2] = fmincon(objective_function_2,x0,A,b,Aeq,beq,lb,ub);   % calling the optimizer
[J(N,1),J(N,2)] = example(x_2);

J(:,1) = J(:,1) + abs(J(1,1));
J(:,2) = J(:,2) + abs(J(N,2));
global a

for i = 1:N
    a = (i-1)/(N-1);
    J(i,1) = a*(J(1,1) + J(N,1));  
    J(i,2) = (1-a)*(J(1,2) + J(N,2));
    
    x0 = [0,0]; % initial point [Nb,p,db]
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [-3 -3];   % Lower bounds
    ub = [3 3];  % Upper bounds

    [xx(i,:),Fval(i,:)] = fmincon(@pareto_objective,x0,A,b,Aeq,beq,lb,ub);   % calling the optimizer
    [obj1,obj2] = example(xx(i,:));
    j(i,1) = obj1;
    j(i,2) = obj2;
end

j(:,1) = j(:,1) + 3;
j(:,2) = j(:,2) + 3;
figure(1)
plot(J(:,1),J(:,2),'*')
hold on
plot(j(:,1),j(:,2),'*')
hold off

function y = pareto_objective(x)
    [obj1,obj2] = example(x);
    global a
    y = a*obj1 + (1-a)*obj2;     % Multi-objective
end

function y = objective_1(x)
    [obj1,obj2] = example(x);
    a = 1;
    y = a*obj1 + (1-a)*obj2;     % Multi-objective
end

function y = objective_2(x)
    [obj1,obj2] = example(x);
    a = 0;
    y = a*obj1 + (1-a)*obj2;     % Multi-objective
end

function [obj1,obj2] = example(x)
    obj1 = (x(1)^2) + x(2);
    obj2 = (x(2)^2) - x(1);
end