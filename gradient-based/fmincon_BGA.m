% Optimization with Gradient Based method - with multistart
clear; clc;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [25 0.0008 0.0006 0.0002];   % Lower bounds
ub = [2500 0.0020 0.0020 0.0020];  % Upper bounds

ConstraintFunction = @Constraints_gradient;

options = optimoptions('fmincon','OutputFcn',{@myplotfval},'Display','iter','Algorithm','sqp');

objective_function = @fmincon_objective; % objective function

for i = 1:10
    x0(i,:) = lb(1,:) + rand.*(ub(1,:) - lb(1,:));
    [x(i,:),Fval(i,:),exitflag(i,:),output(i,:),lambda(i,:),grad(i,:),hessian(:,:,i)] = fmincon(objective_function,x0(i,:),A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
end

optimal = Fval(1,1);
a = 1;
for i = 1:10
    if Fval(i,1) < optimal
        optimal = Fval(i,1);
        a = i;
    end
end

% Global optimal solution
x_opt(1,:) = x(a,:);
Fval_opt(1,1) = Fval(a,1);
grad_opt(1,:) = grad(a,:);
hessian_opt = reshape(hessian(:,:,a),4,4);
eigen_values = eig(hessian_opt);