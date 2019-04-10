%% Optimize params for intensity
clc; clear;
minval = Inf;

lb = [1 1 1 1 1 1 1]/100;
ub = [20 20 20 20 20 20 20];
plb = [3 3 3 3 3 3 3];
pub = [15 15 15 15 15 15 15];

for k = 1:10
    fprintf('Trial # %d\n',k)
    startPoint = datasample(1:20,7);
    %[X, fval] = fminsearch(@costfunction,startPoint);
    [X, fval] = bads(@cost_intCV,startPoint,lb,ub,plb,pub);

    if fval <= minval
        minval = fval;
        params = X;
        sp = startPoint;
    end
end
fprintf('Best values: \n')
[minval params]

% params = [13.7329   10.6318    5.5290    7.6326    4.9132    5.7384    9.2218];
% params = [11.9471    9.9357    3.4262   14.5358    3.1582   13.9762    9.6390];
% params = [14.6602   13.1367   10.7227   14.1914   12.6211    6.3867   10.4883];
% params = [7.2395    5.5119    5.2874    4.7782    4.5480    3.2301   12.1112];
% params = [13.8281    8.6718    4.3592   14.8595   13.2658   13.8280    8.0156];
% params = [11.9774    9.0755   13.2172    8.4293    5.5344    3.7020    3.9625];
% params = [11.8641    9.6456    2.4884    7.5608   10.3759   15.6141   18.7859];
% params = [8.3019   10.7301   12.1513    6.9988   11.9909    9.7962   13.5834];



%% Optimize params for frequency
clc; clear;
minval = Inf;

lb = [1 1 1 1 1 1 1]/100;
ub = [20 20 20 20 20 20 20];
plb = [3 3 3 3 3 3 3];
pub = [15 15 15 15 15 15 15];

for k = 1:10
    fprintf('Trial # %d\n',k)
    startPoint = datasample(1:20,7);
    [X, fval] = bads(@cost_freqCV,startPoint,lb,ub,plb,pub);

    if fval <= minval
        minval = fval;
        params = X;
        sp = startPoint;
    end
end
fprintf('Best values: \n')
[minval params]

% params = [3.4571    5.5904    8.8012    9.6475    3.7618    7.5580    9.6285];
% params = [9.3295    8.1574    7.1621    6.8887    7.4545    7.8117    6.2389];
% params = [5.7331    5.3509    3.6649    9.6311    9.6243    5.7321    9.2928];
% params = [9.3369    8.6926    6.1891    9.0073    2.4794    9.5180    5.4934];
% params = [4.9980    6.9980    4.9980    7.0020    9.9980    9.0020    9.0000];
% params = [4.6417    8.7872    7.8592    7.1048    8.0678    6.6113    7.8465];
% params = [8.0312    9.7739    7.4463   10.5516    6.6093    6.1278   10.2475];
% params = [6.9544   10.4866    3.6779    6.8008    5.4087    5.3772    8.3152];



