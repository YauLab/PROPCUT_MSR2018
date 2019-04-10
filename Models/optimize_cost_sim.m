%% Optimize params for intensity
clc; clear;
minval = Inf;

lb = [1 1 1 1 -100]/100;
ub = [20 20 20 20 1];
plb = [1 1 1 1 -0.1];
pub = [10 10 10 10 0.1];

for k = 1:100
    fprintf('Trial # %d\n',k)
    startPoint = [datasample(1:20,4) rand];
    %[X, fval] = fminsearch(@costfunction,startPoint);
    [X, fval] = bads(@cost_int_sim,startPoint,lb,ub,plb,pub);

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

% params = [7.8127    0.2859    3.9028    4.6204]; %with full
% load('intFullData.mat');

params = [2.7332    4.8790    0.0628    5.0219   -0.0338]; %with 1st half
load('int2ndHalfData.mat')
out.distractor = repmat([0.2 0.4 0.7],1,24)';

% params = [1.1040    1.2331    0.1437    6.8959    0.0627]; %with 2nd half
% load('int1stHalfData.mat')
% out.distractor = repmat([0.2 0.4 0.7],1,24)';

IT = 0.4;
ICv = [0.05 0.20 0.32 0.40 0.48 0.60 0.75];
jndD1 = unique(out.jndD1,'stable');
jndD2 = unique(out.jndD2,'stable');

% Parameters
gamma = params(4);
gc = params(3);
k = params(2);
a = params(1);

for it = 1:10

for ii = 1:72
    ID = out.distractor(ii);
    
    Ihat = (IT + k * IT.^a).^gamma + ( gc * (ID + k * ID.^a) ).^gamma;
    %jndhat = (out.jndD1(ii) + k * out.jndD1(ii).^a).^gamma + ( gc * (out.jndD2(ii) + k * out.jndD2(ii).^a) ).^gamma;

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
        p(ii,jj) = sum( (comparison-target) > 0) / n;
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
% sort(mean(r2,2),'descend')
meanr2(it,:) = mean(r2,2)';
mean(mean(r2,2));
end

% r2 for each half. Use these values to generate figure for variance
% explained for each subject and within subject CV.
r2 = [0.9155 0.9467 0.9238 0.9248 0.9109 0.9225 0.9160 0.9407;
      0.9255 0.9414 0.9355 0.9507 0.9092 0.9165 0.9162 0.9537];
meanr2 = mean(r2);
stdr2 = std(r2)
sort(meanr2)
mean(meanr2)



%% Optimize params for intensity
clc; clear;
minval = Inf;

lb = [1 1 1 1 -100]/100;
ub = [20 20 20 20 1];
plb = [1 1 1 1 -0.1];
pub = [10 10 10 10 0.1];

for k = 1:100
    fprintf('Trial # %d\n',k)
    startPoint = [datasample(1:20,4) rand];
    %[X, fval] = fminsearch(@costfunction,startPoint);
    [X, fval] = bads(@cost_freq_sim,startPoint,lb,ub,plb,pub);

    if fval <= minval
        minval = fval;
        params = X;
        sp = startPoint;
    end
    [minval params]
end
fprintf('Best values: \n')
[minval params]
