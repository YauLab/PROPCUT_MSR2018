function cost = cost_freq(params)

load('freq1stHalfData.mat');
out.distractor = repmat([100 200 300],1,24)';

FT = 200;
FCv = [100 140 180 200 220 260 300];
jndD1 = unique(out.jndD1,'stable');

% Parameters
alpha = params(1);
beta = params(2);
gamma = params(3);
ge = params(4);
gc = params(5);
k = params(6);
a = params(7);

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
    
    n = 100;
    for jj = 1:7
        FC = FCv(jj);
        compWeights = normpdf(FCv,FC,jndD1(subj));
        compWeights = compWeights / sum(compWeights);
        comparison = datasample(FCv,n,'Replace',true,'Weights',compWeights);
        target = datasample(FCv,n,'Replace',true,'Weights',targWeights);
        p(ii,jj) = sum(comparison>target) / n;
    end
end
cost = sum((out.choiceprob(:)-p(:)).^2);
