function [state,options,optchanged] = gahistory(options,state,flag)
     global d1 it
    optchanged = false;
    switch flag
    %case 'init'
    %    h1 = figure;
    %    history(:,:,1) = state.Population;
    %    assignin('base','gapopulationhistory',history);
    case 'iter'
        if state.Generation == 1
             d1 = [];
             it = [];
        end
        y = state.Best;
        j = state.Generation;
        it(j) = j;
        d1 = y;
    case 'done'
%         plot(it,d1,'*-')
%         title('Optimizer History plot')
%         xlabel('Generation number')
%         ylabel('fitness value')
%         legend('Run 1','Run 2','Run 3','Run 4','Run 5','Run 6','Run 7','Run 8','Run 9','Run 10')       
%         grid on
%         hold on

   end