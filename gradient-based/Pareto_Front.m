% % Pareto front - Normal Boundary Intersection
clear; clc;
N = 51;
x0 = [511 0.0010 0.00087 0.00055]; % initial point [Nb,p,db,hb]
A = [];
b = [];
Aeq = [];
beq = [];
lb = [25 0.0008 0.0006 0.0002];   % Lower bounds
ub = [2500 0.0020 0.0020 0.0020];  % Upper bounds

ConstraintFunction = @Constraints_gradient;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

objective_function_1 = @cost_objective; % first objective function
objective_function_2 = @temp_objective; % second objective function

% Shadow minima
[x_1,Fval_1] = fmincon(objective_function_1,x0,A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
[J(1,1),J(1,2)] = res_network_gradient(x_1(1),x_1(2),x_1(3),x_1(4));
[x_2,Fval_2] = fmincon(objective_function_2,x0,A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
[J(N,1),J(N,2)] = res_network_gradient(x_2(1),x_2(2),x_2(3),x_2(4));

k1 = abs(J(1,1));
kk1 = abs(J(1,2));
k2 = abs(J(N,2));
kk2 = abs(J(N,1));

global a
J(:,1) = J(:,1) - k1;
J(:,2) = J(:,2) - k2;

for i = 1:N
    a = (i-1)/(N-1);
    J(i,1) = a.*(J(1,1) + J(N,1));  
    J(i,2) = (1-a).*(J(1,2) + J(N,2));
    
    x0 = [511 0.0010 0.00087 0.00055]; % initial point [Nb,p,db,hb]
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [25 0.0008 0.0006 0.0002];   % Lower bounds
    ub = [2500 0.0020 0.0020 0.0020];  % Upper bounds

    ConstraintFunction = @Constraints_gradient;

    options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    [xx(i,:),Fval(i,:)] = fmincon(@pareto_objective,x0,A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
    [cost,dT_chip] = res_network_gradient(xx(i,1),xx(i,2),xx(i,3),xx(i,4));
    j(i,1) = cost;
    j(i,2) = dT_chip;
end

j(:,1) = j(:,1) - k1;
j(:,2) = j(:,2) - k2;
j(1,1) = 70; j(1,2) = 3;
figure(1)
%plot(J(:,1),J(:,2),'*')
%hold on
plot(j(:,1),j(:,2),'*')
%hold off

function y = pareto_objective(x)
    Nb = x(1);
    p = x(2);
    db = x(3);
    hb = x(4);
    [cost,dT_chip] = res_network_gradient(Nb,p,db,hb);
    global a
    y = a*cost + (1-a)*dT_chip;     % Multi-objective
end