%% Generate INFO 2 iteration sequence table
clear input;
clc;
% Add data, change as needed
idx = 1:length(info2.seq');
input.data = [idx(:) info2.seq']; % make sure to have this correctly formatted
input.tablePositioning = 'h';
input.tableColLabels = {'Iteration', '$x_1$','$x_2$','$x_3$','$x_4$','$x_5$'};
input.dataFormat = {'%d',1,'%.4f',5};
input.tableColumnAlignment = 'l';
input.tableCaption = 'Iteration sequence for BFGS approximation';
input.tableLabel = 'ex2:bfgs'; % prepends "table" -> table:MyTableLabel
% Generate tex output
texOut = latexTable(input);

%% Generate INFO 3 iteration sequence table
clear input;
clc;
% Add data, change as needed
idx = 1:length(info3.seq');
input.data = [idx(:) info3.seq']; % make sure to have this correctly formatted
input.tablePositioning = 'h';
input.tableColLabels = {'Iteration', '$x_1$','$x_2$','$x_3$','$x_4$','$x_5$'};
input.dataFormat = {'%d',1,'%.4f',5};
input.tableColumnAlignment = 'l';
input.tableCaption = 'Iteration sequence for BFGS approximation and line search';
input.tableLabel = 'ex2:bfgsls'; % prepends "table" -> table:MyTableLabel
% Generate tex output
texOut = latexTable(input);