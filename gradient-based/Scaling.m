% Assignment 4 B

% Scaling
clear
clc

% Calculation of Hessian

A = [];
b = [];
Aeq = [];
beq = [];
lb = [25 0.0008 0.0006 0.0002];   % Lower bounds
ub = [2500 0.0020 0.0020 0.0020];  % Upper bounds

ConstraintFunction = @Constraints_gradient;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

objective_function = @fmincon_objective; % objective function

x0(1,:) = [511 0.0010 0.00087 0.00055];
tic
[x(1,:),Fval(1,:),exitflag(1,:),output(1,:),lambda(1,:),grad(1,:),hessian(:,:,1)] = fmincon(objective_function,x0(1,:),A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
time = toc;

% Global optimal solution
x_opt(1,:) = x(1,:);
Fval_opt(1,1) = Fval(1,:);
%grad_opt(1,:) = grad(1,:);
hessian_opt = reshape(hessian(:,:,1),4,4);
eigen_values = eig(hessian_opt);

dy1 = diag(hessian_opt);
x1 = x_opt;
dx1 = [0.0001 0.00001 0.00001 0.00001];
disp('Diagonal values of Hessian are');disp(dy1)

% Calculating the scaling factor for each design variable

L(1) = 1500;
L(2) = 0.0004*(1/dy1(2))^(0.02);
L(3) = 0.0009*(1/dy1(3))^(0.03);
L(4) = 0.0001*(1/dy1(4))^(0.03);

disp('Scaling factor for variables are');disp(L)

% Updated variables and function

x2 = x1./L;
y2 = res_network_gradient(x2(1),x2(2),x2(3),x2(4));
disp('New Fvalue');disp(y2)
disp('New scaled x* =');disp(x2)

% Re-calculate the Hessian after scaling

dx2 = dx1./L;

xx2(1,:,:) = [L(1)*(x2(1)-dx2(1))   L(2)*x2(2)          L(3)*x2(3)          L(4)*x2(4)
              L(1)*x2(1)            L(2)*x2(2)          L(3)*x2(3)          L(4)*x2(4)
              L(1)*(x2(1)+dx2(1))   L(2)*x2(2)          L(3)*x2(3)          L(4)*x2(4)];
xx2(2,:,:) = [L(1)*x2(1)            L(2)*(x2(2)-dx2(2)) L(3)*x2(3)          L(4)*x2(4)
              L(1)*x2(1)            L(2)*x2(2)          L(3)*x2(3)          L(4)*x2(4)
              L(1)*x2(1)            L(2)*(x2(2)+dx2(2)) L(3)*x2(3)          L(4)*x2(4)];
xx2(3,:,:) = [L(1)*x2(1)            L(2)*x2(2)          L(3)*(x2(3)-dx2(3)) L(4)*x2(4)
              L(1)*x2(1)            L(2)*x2(2)          L(3)*x2(3)          L(4)*x2(4)
              L(1)*x2(1)            L(2)*x2(2)          L(3)*(x2(3)+dx2(3)) L(4)*x2(4)];
xx2(4,:,:) = [L(1)*x2(1)            L(2)*x2(2)          L(3)*x2(3)          L(4)*(x2(4)-dx2(4))
              L(1)*x2(1)            L(2)*x2(2)          L(3)*x2(3)          L(4)*x2(4)
              L(1)*x2(1)            L(2)*x2(2)          L(3)*x2(3)          L(4)*(x2(4)+dx2(4))];
%disp(xx2(1,1,:))
%disp(x2.*L)
for i = 1:4
    dy2(i) = (res_network_gradient(xx2(i,1,1),xx2(i,1,2),xx2(i,1,3),xx2(i,1,4)) - 2*res_network_gradient(xx2(i,2,1),xx2(i,2,2),xx2(i,2,3),xx2(i,2,4)) + res_network_gradient(xx2(i,3,1),xx2(i,3,2),xx2(i,3,3),xx2(i,3,4)))/(dx2(i)^2);
end

disp('Diagonal values of Hessian after scaling are - ');disp(dy2)

% Optimization after scaling
A = [];
b = [];
Aeq = [];
beq = [];
lb = [25 0.0008 0.0006 0.0002]./L;   % Lower bounds
ub = [2500 0.0020 0.0020 0.0020]./L;  % Upper bounds

ConstraintFunction = @Constraints_gradient;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

objective_function = @fmincon_scaled_objective; % objective function
x0(1,:) = [511 0.0010 0.00087 0.00055]./L;
tic
[x(1,:),Fval(1,:),exitflag(1,:),output(1,:),lambda(1,:),grad(1,:),hessian(:,:,1)] = fmincon(objective_function,x0(1,:),A,b,Aeq,beq,lb,ub,ConstraintFunction,options);   % calling the optimizer
time_scaled = toc;
x_scaled_opt(1,:) = L.*x(1,:);