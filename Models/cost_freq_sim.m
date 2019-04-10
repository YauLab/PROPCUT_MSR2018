function cost = cost_freq_sim(params)

load('freq1stHalfData.mat');
out.distractor = repmat([100 200 300],1,24)';

FT = 200/1000;
FCv = [100 140 180 200 220 260 300]/1000;
jndD1 = unique(out.jndD1,'stable');

% Parameters
gamma = params(4);
gc = params(3);
k = params(2);
a = params(1);

for ii = 1:72
    FD = out.distractor(ii)/1000;
    
    Fhat = (FT + k * FT.^a).^gamma + ( gc * (FD + k * FD.^a) ).^gamma;
    
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
        p(ii,jj) = sum(comparison-target+params(5)) / n;
    end
end
cost = sum((out.choiceprob(:)-p(:)).^2);
