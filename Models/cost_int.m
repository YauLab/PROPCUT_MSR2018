function cost = cost_int(params)

load('int1stHalfData.mat')
out.distractor = repmat([0.2 0.4 0.7],1,24)';

IT = 0.4;
ICv = [0.05 0.20 0.32 0.40 0.48 0.60 0.75];
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
cost = sum((out.cp(:)-p(:)).^2);
