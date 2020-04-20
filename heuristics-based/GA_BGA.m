%GA Application
clear; clc;
ObjectiveFunction = @ga_objective;
nvars = 4;    % Number of variables
lb = [25 0.0008 0.0006 0.0002];   % Lower bounds
ub = [2500 0.0020 0.0020 0.0020];  % Upper bounds
ConstraintFunction = @Constraints_heuristic;
%PopulationSize = 500;
IntCon = 1;
%options = optimoptions(@ga,'PlotFcn',@gahistory,'Display','iter','PopulationSize',500,'UseParallel',true);
    options = optimoptions(@ga,'PlotFcn',{@gahistory},'Display','iter',...
                           'PopulationSize',100,...%'MutationFcn',{@mutationadaptfeasible},...
                           'InitialPenalty',1100,...
                           'EliteCount',1,...
                           'CrossoverFraction',0.9,...
                           'MaxGenerations',100,...
                           'UseParallel',true);
                       
    options.FunctionTolerance=1e-16;
    options.ConstraintTolerance=1e-8;
%for i = 1:10
    [x(1,:),fval(1,:)] = ga(ObjectiveFunction,nvars,[],[],[],[],lb,ub,ConstraintFunction,IntCon,options);
    global d1 it
    plot(it,d1,'*-')
    title('Optimizer History plot')
    xlabel('Generation number')
    ylabel('fitness value')
    legend('Run 1','Run 2','Run 3','Run 4','Run 5','Run 6','Run 7','Run 8','Run 9','Run 10')       
    grid on
%end


