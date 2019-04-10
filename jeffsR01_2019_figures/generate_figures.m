clc; clear;
% % Set some cool plot/graphics environments!
% c =[ 0 0 1;...
%      1 0 1;...
%      1 0 0];
% set(0,'defaultAxesColorOrder',c);
% set(0,'defaultlinelinewidth',2)
% set(0,'DefaultAxesFontSize',16)

% Read data as a structure
load('freqFullData.mat')

% Values used in the experiment
standard = 200;
% comparison frequency
fc = repmat([100 140 180 200 220 260 300], 9, 1);
% hand positions
prop = repmat([0*ones(3,1); 0.5*ones(3,1); ones(3,1)], 1, 7);
% distractor frequency
fd = repmat( repmat([100;200;300],3,7), 1, 1);
% number of repeats
n = repmat([20 22 24 24 24 22 20], 9, 1);

% Initialization for parameter search
c_init_si = [
    0.4524      -0.8804      0.49372      0.95277        18097
    0.4695       -0.763      0.50442       1.0561        18446
    0.4793      -1.2629      0.53364      0.56515        13806
    0.5         -1.625       0.52904      0.52594        14026
    0.4539      -0.7704      0.48656       1.0659        22936
    0.4541      -1.0743      0.52087      0.69247        17071
    0.4671      -1.3744      0.45694      0.59343        18429
    0.4693      -0.7823      0.55619       1.0338        17290];
%
% 2D frequency- and location-based modulation function
afs = @(x1,x2,x3,dp,df) x1 * exp(-( dp/x2 + abs(df).^2/x3 ));
afsn = @(x1,x2,dp,df) exp(-( dp/x1 + abs(df).^2/x2 ));
% 1D frequency-based modulation function
afsf = @(x1,x2,df) x1 * exp(-( abs(df).^2/x2 ));
afsfn = @(x1,df) exp(-( abs(df).^2/x1 ));
% 1D location-based modulation function
afsp = @(x1,x2,dp) x1 * exp(-( dp/x2 ));
afspn = @(x1,dp) exp(-( dp/x1 ));

for subj = 1:8
    %
    % Group CV
    fc = repmat([100 140 180 200 220 260 300],9*7,1);
    prop = repmat([zeros(3,1); 0.5*ones(3,1); ones(3,1)],1*7,7);
    fd = repmat( repmat([100;200;300],3,7), 7,1);
    n = repmat([20 22 24 24 24 22 20],9*7,1);
    
    jndD1 = repmat(out.jndD1,1,7); jndD1(9*subj-8:9*subj,:) = [];
    jndD2 = repmat(out.jndD2,1,7); jndD2(9*subj-8:9*subj,:) = [];
    subj_cp = out.choiceprob; subj_cp(9*subj-8:9*subj,:) = [];
    %}
    
    %
    % Noise ceiling
    nspse = out.pse; nspse(9*subj-8:9*subj)=[];
    nsjnd = out.jnd; nsjnd(9*subj-8:9*subj)=[];
    mycp = reshape(subj_cp',7,9,7);
    mcp_ns(9*subj-8:9*subj,:) = mean(mycp,3)';
    pse_ns(9*subj-8:9*subj,1) = mean(reshape(nspse,9,7),2);
    jnd_ns(9*subj-8:9*subj,1) = mean(reshape(nsjnd,9,7),2);
    %}

    
    % Count of subject choice probability
    X = subj_cp .* n;

    
    % F+P, with normalization
    % 1Df, 1Dp - ra
    wT = @(c,prop) 1 ./ ( c(4)+afspn(c(3),prop) );
    wD = @(c,prop) afspn(c(3),prop) ./ ( c(4)+afspn(c(3),prop) );
    %wT = @(c,df,jndD1,jndD2) 1 ./ ( (sqrt(jndD1.^2+jndD2.^2)/(c(1)*jndD1))+afsfn(c(3),df) );
    %wD = @(c,df,jndD1,jndD2) afsfn(c(3),df) ./ ( (sqrt(jndD1.^2+jndD2.^2)/(c(1)*jndD1))+afsfn(c(3),df) );
    jndd1 = @(c,jndD1) jndD1 + c(5);
    jndd2 = @(c,jndD2) jndD2 + c(6);
    f1Df1Dpra = @(c,fc,prop,standard,fd,jndD1,jndD2)...
        -sum(sum(...
        X.*log(0.5 * ( 1 + erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndd1(c,jndD1)/c(1)).^2 + (wD(c,prop).*jndd2(c,jndD2)./afsf(c(1),c(2),(standard-fd))).^2 ) )...
        ) ) ) + ...
        (n-X).*log(0.5 * ( 1 - erf(...
        (fc - wT(c,prop)*standard - wD(c,prop).*fd ) ./...
        sqrt( 2 * ( (wT(c,prop).*jndd1(c,jndD1)/c(1)).^2 + (wD(c,prop).*jndd2(c,jndD2)./afsf(c(1),c(2),(standard-fd))).^2 ) )...
        ) ) )...
        ));
    jndx1 = unique(out.jndD1); 
    jndx2 = unique(out.jndD2);
    fun = @(c)f1Df1Dpra(c,fc,prop,standard,fd,jndD1,jndD2);
    c1Df1Dpra(subj,:) = fminsearch( fun, [c_init_si(subj,[3 5 4]) 1 0 0], optimset('MaxFunEvals',100) );
    ll_DfDpra(subj) = -f1Df1Dpra(c1Df1Dpra(subj,:),fc,prop,standard,fd,jndD1,jndD2);

end
%%
%
clc;

% Compute model probabilities based on estimated c parameters
prop = 0.5;
row = 0;
for subj = 1:8
    pseD1 = out.pseD1(9*subj - 8);
    jndD1 = out.jndD1(9*subj - 8);
    jndD2 = out.jndD2(9*subj - 8);
    for fd = [100 200 300]
        row = row + 1;
        wT = 1 ./ ( c1Df1Dpra(subj,4)+afspn(c1Df1Dpra(subj,3),prop) );
        wD = afspn(c1Df1Dpra(subj,3),prop) ./ ( c1Df1Dpra(subj,4)+afspn(c1Df1Dpra(subj,3),prop) );
        pse_DfDpra(row) = wT*standard + wD.*fd;
        jnd_DfDpra(row) = sqrt( (wT*jndD1/c1Df1Dpra(subj,1)).^2 + (wD*jndD2./afsf(c1Df1Dpra(subj,1),c1Df1Dpra(subj,2),(standard-fd))).^2 );
    end
end

%%
clc;
cmap = parula;
cmap(1,:) = [0.7 0.7 0.7];

f_standard = linspace(100,340,5);
fd = [100, 220, 340];


labels = {'L100','L130','L160','L190','L220','L250','L280','L310','L340'...
    'R100','R130','R160','R190','R220','R250','R280','R310','R340'...
    'L100-R100','L220-R100','L340-R100','L100-R160','L220-R160','L340-R160'...
    'L100-R220','L220-R220','L340-R220','L100-R280','L220-R280','L340-R280'...
    'L100-R340','L220-R340','L340-R340'};


prop = 0.5;
c = mean(c1Df1Dpra);

wT = 1 ./ ( c(4)+afspn(c(3),prop) );
wD = afspn(c(3),prop) ./ ( c(4)+afspn(c(3),prop) );

com = [];
for standard = f_standard
    com(end+1:end+3) = wT*standard + wD.*fd;
end

avg = [];
for standard = f_standard
    avg(end+1:end+3) = 0.5*standard + 0.5*fd;
end

% Unimanual
ur = linspace(100,340,9);
ul = linspace(100,340,9);
% Bimanual, right
br = repmat(f_standard',1,3)'; br = br(:)';
% Bimanual, left
bl = repmat(fd, 1, 5);
% Bimanual, combined
bc = com;

% Bimanual, average
ba = avg;

colors = [0.7, 0.7, 0.7];

% Left
figure(1)
x = [ul, ur, bl];
d = pdist(x');
d = squareform(d);
for i = 1:33
    for j = 10:18
        if i ~= j
            d(i,j) = -10;
        end
    end
end
for i = 10:18
    for j = 1:33
        if i ~= j
            d(i,j) = -10;
        end
    end
end
for i = 1:33
    for j = 1:33
        if i < j
            d(i,j) = nan;
        end
    end
end
[nr,nc] = size(d);
pcolor([d nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gca,'XTick',1.5:33.5,'XTickLabel',labels)
set(gca,'YTick',1.5:33.5,'YTickLabel',labels)
xtickangle(90)
axis image
set(gca,'TickDir','out')
set(gca,'color',[0.2 0.2 0.2])
colormap(cmap)
grid('off')



% Right
figure(2)
x = [ul, ur, br];
d = pdist(x');
d = squareform(d);
for i = 1:33
    for j = 1:9
        if i ~= j
            d(i,j) = -10;
        end
    end
end
for i = 1:9
    for j = 1:33
        if i ~= j
            d(i,j) = -10;
        end
    end
end
for i = 1:33
    for j = 1:33
        if i < j
            d(i,j) = nan;
        end
    end
end
[nr,nc] = size(d);
pcolor([d nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gca,'XTick',1.5:33.5,'XTickLabel',labels)
set(gca,'YTick',1.5:33.5,'YTickLabel',labels)
xtickangle(90)
axis image
set(gca,'TickDir','out')
set(gca,'color',[0.2 0.2 0.2])
colormap(cmap)


% Combined, propcut model
figure(3)
x = [ul, ur, bc];
d = pdist(x');
d = squareform(d);
for i = 1:33
    for j = 1:33
        if i < j
            d(i,j) = nan;
        end
    end
end
[nr,nc] = size(d);
pcolor([d nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gca,'XTick',1.5:33.5,'XTickLabel',labels)
set(gca,'YTick',1.5:33.5,'YTickLabel',labels)
xtickangle(90)
axis image
set(gca,'TickDir','out')
set(gca,'color',[0.2 0.2 0.2])


% Combined, average
figure(4)
x = [ul, ur, ba];
d = pdist(x');
d = squareform(d);
for i = 1:33
    for j = 1:33
        if i < j
            d(i,j) = nan;
        end
    end
end
[nr,nc] = size(d);
pcolor([d nan(nr,1); nan(1,nc+1)]);
shading flat;
set(gca, 'ydir', 'reverse');
set(gca,'XTick',1.5:33.5,'XTickLabel',labels)
set(gca,'YTick',1.5:33.5,'YTickLabel',labels)
xtickangle(90)
axis image
set(gca,'TickDir','out')
set(gca,'color',[0.2 0.2 0.2])
