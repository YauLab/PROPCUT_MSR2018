function cost = cost_int_sim(params)

load('int2ndHalfData.mat')
out.distractor = repmat([0.2 0.4 0.7],1,24)';

IT = 0.4;
ICv = [0.05 0.20 0.32 0.40 0.48 0.60 0.75];
jndD1 = unique(out.jndD1,'stable');

% Parameters
gamma = params(4);
gc = params(3);
k = params(2);
a = params(1);

for ii = 1:72
    ID = out.distractor(ii);

    Ihat = (IT + k * IT.^a).^gamma + ( gc * (ID + k * ID.^a) ).^gamma;

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
cost = sum((out.cp(:)-p(:)).^2);
