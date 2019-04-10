%% Optimize params for frequency
clc; clear;
minval = Inf;

lb = [1 1 1 1 1 1 1]/100;
ub = [20 20 20 20 20 20 20];
plb = [1 1 1 1 1 1 1];
pub = [19 10 10 10 10 10 10];
% 18.1586   15.9945    6.9929    7.9871    2.0013    7.0038    5.0125    9.2077
for k = 1:10
    fprintf('Trial # %d\n',k)
    startPoint = datasample(1:20,7);
    %[X, fval] = fminsearch(@costfunction_freq,startPoint);
    [X, fval] = bads(@cost_freq,startPoint,lb,ub,plb,pub);

    if fval <= minval
        minval = fval;
        params = X;
        sp = startPoint;
    end
end
fprintf('Best values: \n')
[minval params]

%% Analyze frequency data based on optimized params
clc; clear;
% Parameters

% params = [9.8062    4.4708    2.2571    2.1577    4.8053    1.7545    9.6901]; %with full data
% load('freqFullData.mat');
% out.distractor = repmat([100 200 300],1,24)';

% Params estimated with 1st half and test the model with 2nd half
params = [1.0726    1.9665    3.4577    3.9341    4.8060    5.1271    4.5186];
load('freq2ndHalfData.mat');
out.distractor = repmat([100 200 300],1,24)';

% params = [1.0482    2.2507    6.9041    0.9849    2.4815    3.0852    1.9178]; %with 2nd half
% load('freq1stHalfData.mat');
% out.distractor = repmat([100 200 300],1,24)';



alpha = params(1);
beta = params(2);
gamma = params(3);
ge = params(4);
gc = params(5);
k = params(6);
a = params(7);

FT = 200;
FCv = [100 140 180 200 220 260 300];
jndD1 = unique(out.jndD1);

for ii = 1:72
    FD = out.distractor(ii);
    
    % Energies
    ED = (FD/gc).^gamma;
    ET = (FT/gc).^gamma;
    ETstar = (FT/ge).^gamma;
    EDstar = (FD/ge).^gamma;
    % Target
    FTnum = 1 + EDstar / (1 + beta^gamma * ET);
    FTden = 1 + ED / (1 + alpha^gamma * ET);
    FTMN = FT + k * FT.^a;
    FThat = FTnum ./ FTden + FTMN;
    % Distractor
    FDnum = 1 + ETstar / (1 + beta^gamma * ED);
    FDden = 1 + ET / (1 + alpha^gamma * ED);
    FDMN = FD + k * FD.^a;
    FDhat = FDnum ./ FDden + FDMN;
    % Target + Distractor
    Fhat = FThat + FDhat;
    
    subj = ceil(ii/9);
    
    targWeights = normpdf(FCv,Fhat,jndD1(subj))+eps;
    targWeights = targWeights / sum(targWeights);
    
    n = 10000;
    for jj = 1:7
        FC = FCv(jj);
        compWeights = normpdf(FCv,FC,jndD1(subj));
        compWeights = compWeights / sum(compWeights);
        comparison = datasample(FCv,n,'Replace',true,'Weights',compWeights);
        target = datasample(FCv,n,'Replace',true,'Weights',targWeights);
        p(ii,jj) = sum(comparison>target) / n;
    end
end

for subj = 1:8
    startRow = 9 * subj - 8;
    endRow = 9 * subj;
    pmes = out.choiceprob(startRow:endRow,:);
    pest = p(startRow:endRow,:);
    for cond = 1:9
        r2(subj,cond) = corr(pmes(cond,:)',pest(cond,:)') ^ 2;
    end
end
sort(mean(r2,2),'descend')
mean(mean(r2,2))

%% Optimize params for intensity
clc; clear;
minval = Inf;

lb = [1 1 1 1 1 1 1]/100;
ub = [20 20 20 20 20 20 20];
plb = [1 1 1 1 1 1 1];
pub = [15 15 15 15 15 15 15];
% 21.6308    8.3832    9.4865    9.0551   10.7455    6.0349   14.1930    9.8039
for k = 1:10
    fprintf('Trial # %d\n',k)
    startPoint = datasample(1:20,7);
    %[X, fval] = fminsearch(@costfunction,startPoint);
    [X, fval] = bads(@cost_int,startPoint,lb,ub,plb,pub);

    if fval <= minval
        minval = fval;
        params = X;
        sp = startPoint;
    end
end
fprintf('Best values: \n')
[minval params]


%% Analyze intensity data based on optimized parameters
clear;
% Parameters

% params = [8.3832    9.4865    9.0551   10.7455    6.0349   14.1930    9.8039]; %with full
% load('intFullData.mat');

params = [7.9468    2.7433    7.9672    2.4761   12.1628    1.1994    4.3477]; %with 1st half
load('int2ndHalfData.mat')
out.distractor = repmat([0.2 0.4 0.7],1,24)';

params = [1.7364    2.2288    0.9953    4.9690   18.9238    6.9072    0.6259]; %with 2nd half
load('int1stHalfData.mat')
out.distractor = repmat([0.2 0.4 0.7],1,24)';

alpha = params(1);
beta = params(2);
gamma = params(3);
ge = params(4);
gc = params(5);
k = params(6);
a = params(7);

IT = 0.4;
ICv = [0.05 0.20 0.32 0.40 0.48 0.60 0.75];
jndD1 = [0.08 0.02 0.06 0.08 0.06 0.12 0.09 0.07];

for ii = 1:72
    ID = out.distractor(ii);
    
    % Energies
    ED = (ID/gc).^gamma;
    ET = (IT/gc).^gamma;
    ETstar = (IT/ge).^gamma;
    EDstar = (ID/ge).^gamma;
    % Target
    ITnum = 1 + EDstar / (1 + beta^gamma * ET);
    ITden = 1 + ED / (1 + alpha^gamma * ET);
    ITMN = IT + k * IT.^a;
    IThat = ITnum ./ ITden + ITMN;
    % Distractor
    IDnum = 1 + ETstar / (1 + beta^gamma * ED);
    IDden = 1 + ET / (1 + alpha^gamma * ED);
    IDMN = ID + k * ID.^a;
    IDhat = IDnum ./ IDden + IDMN;
    % Target + Distractor
    Ihat = IThat + IDhat;
    
    
    
    subj = ceil(ii/9);
    
    targWeights = normpdf(ICv,Ihat,jndD1(subj))+eps;
    targWeights = targWeights / sum(targWeights);
    
    n = 100;
    for jj = 1:7
        IC = ICv(jj);
        compWeights = normpdf(ICv,IC,jndD1(subj));
        compWeights = compWeights / sum(compWeights);
        comparison = datasample(ICv,n,'Replace',true,'Weights',compWeights);
        target = datasample(ICv,n,'Replace',true,'Weights',targWeights);
        p(ii,jj) = sum(comparison>target) / n;
    end
end

for subj = 1:8
    startRow = 9 * subj - 8;
    endRow = 9 * subj;
    pmes = out.cp(startRow:endRow,:);
    pest = p(startRow:endRow,:);
    for cond = 1:9
        r2(subj,cond) = corr(pmes(cond,:)',pest(cond,:)') ^ 2;
    end
%     r(subj) = corr(pmes(:),pest(:));
end
sort(mean(r2,2),'descend')
mean(mean(r2,2))
