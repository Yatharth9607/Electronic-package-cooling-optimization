function stop = myplotfval(fval,optimValues,state)
     persistent d1 it
     if ~nargin % example reset mechanism
         d1=[];
     end
     if strcmp(state,'iter')
         if optimValues.iteration == 0
             d1 = [];
             it = [];
         end
         y = optimValues.fval;
         j = optimValues.iteration + 1;
         it(j) = j;
         d1(j) = y;
     end
     if strcmp(state,'done')
         plot(it,d1,'*-')
         title('Optimizer History')
         xlabel('Number of Iterations')
         ylabel('Function value')
         legend('Run 1','Run 2','Run 3','Run 4','Run 5','Run 6','Run 7','Run 8','Run 9','Run 10')       
         grid on
         hold on
        %filename = 'Golinski.xlsx';
        %xlswrite(filename,d1)
     end
     stop = 0;
     

   end