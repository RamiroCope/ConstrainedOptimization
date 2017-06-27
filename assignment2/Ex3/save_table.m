function save_table(table_name,x0,iter1,iter2,iter3,fcalls1,fcalls2,fcalls3,T1,T2,T3)

T1 = 1000*T1; 
T2 = 1000*T2;
T3 = 1000*T3


??file = fopen(['../table/' table_name],'w');
fprintf(file,'\\begin{tabular}{l | ccc | ccc | ccc} \\hline \\hline
fprintf(file,'& \\multicolumn{3}{c}{$x 0 = (%d,%d)$} & \\multicolumn{3}{c}{$x 0 = (%d,%d)$} & \ fprintf(file,'Method & Iter & Evals & Time & Iter & Evals & Time & Iter & Evals & Time \\\\ \\h
fprintf(file,'Backtracking & %d & %d & %.1f & %d & %d & %.1f & %d & %d & %.1f \\\\ \n',iter1(1) fprintf(file,'Soft line search & %d & %d & %.1f & %d & %d & %.1f & %d & %d & %.1f \\\\ \n',iter
fprintf(file,'\\hline \\hline \n'); fprintf(file,'\\end{tabular} \n'); fclose(file);

