
%{

--- ClusterAnalysis ---
Loads the data from a series of simulations runs, across different seeds,
and combines all the data, normalizes and clusters them using a k-means
clustering. 

This scripts expects a series of subsequent simulation runs, with saveIDs
starting with an 'M' and followed by a number. These numbers first loop
through different seeds, then different motifs. For example, when three seeds 
were used, run M1 used seed 1 for motif I, M2 seed 2 for motif I, M3 seed 3 
for motif I and M4 seed 1 for motif II, etc. The motifs and seeds that
should be included for clustering can be specified in lines 36 and 37.

Note that the k-means clustering depends on random initial conditions;
clustering can therefore produce different cluster orders on
subsequent passes of the code. 

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

clear all

%% parameters
path_settings = './scripts_settings/';
path_data = 'YourPathHere/Analysis/';
path_figures = 'YourPathHere/Figures';

parfiles = dir([path_settings,'*_Settings.m']);
firstSaveID = 1; % the first saveID from the series

% data to use
parincl = [1:20]; %[1,2];
seedincl = 1:10; % seeds to include in the analysis

seeds = 1:10; % all seeds used in the runs (in case seeds need to be skipped)

%% load data and prepare for clustering and plot analysis

for pp = parincl
    fprinft('Loading motif %i\n',pp)
    seednr = seeds(1);
    % determine the run's saveID
    saveID = ['M', num2str(firstSaveID + (pp-1)*length(seeds)+seedincl(1))];
    feval(parfiles(pp).name(1:end-2));
    
    for ss = seedincl
        
        saveID = ['M', num2str(firstSaveID + (pp-1)*length(seeds)+ss)];
        seednr = seeds(ss);
        
        % load
        load([path_datasaveID,'_Analysis.mat']);
        allvars = whos('-file',[path_datasaveID,'_Analysis.mat']);
        
        % compute sum across seeds
        if seednr == seedincl(1)
            for v = 1:length(allvars)
                if iscell(eval([allvars(v).name,';'])) || strcmp(allvars(v).name(1:3), 'fce') 
                    continue;  
                end
                eval([genvarname([allvars(v).name,'_avg']), '=' allvars(v).name,';']);
                eval([genvarname([allvars(v).name,'_size']), '= ~isnan(' allvars(v).name,');']);
            end
        else
            for v = 1:length(allvars)
                if iscell(eval([allvars(v).name,';'])) || strcmp(allvars(v).name(1:3), 'fce')
                    continue;  
                end
                eval([genvarname([allvars(v).name,'_avg']),' = nansum(cat(3,',...
                    genvarname([allvars(v).name,'_avg'])...
                    ',' allvars(v).name,'),3);']);
                eval([genvarname([allvars(v).name,'_size']),' = ',...
                    genvarname([allvars(v).name,'_size'])...
                    '+ ~isnan(' allvars(v).name,');']);
            end
        end
        clearvars(allvars.name)
    end
   
    % average across seeds
    for v = 1:length(allvars)
        if ~exist(genvarname([allvars(v).name, '_avg'])); continue;  end
        eval([genvarname([allvars(v).name,'_avg']),' = ',...
            genvarname([allvars(v).name,'_avg']) './',...
            genvarname([allvars(v).name,'_size']),';']);
        
        % exclude if fewer than 5 seeds contain a value
        % these values will be set to 0 after normalization and before
        % clustering - see line 138
        tmp = eval(genvarname([allvars(v).name,'_size']))<=5;
        eval([genvarname([allvars(v).name,'_avg']),'(tmp==1) = NaN;']);
    end
    
    % prep for clustering
    if pp == parincl(1)
        for v = 1:length(allvars)
            if ~exist(genvarname([allvars(v).name, '_avg'])); continue;  end
            eval([genvarname([allvars(v).name,'_clas']),' = ',...
                genvarname([allvars(v).name,'_avg']), ';']);
        end
    else
        for v = 1:length(allvars)
            if ~exist(genvarname([allvars(v).name, '_avg'])); continue;  end
            if size(eval(genvarname([allvars(v).name,'_clas'])),1) >...
                    size(eval(genvarname([allvars(v).name,'_avg'])),1)
                dum = nan(size(eval(genvarname([allvars(v).name,'_clas'])),1),...
                    size(eval(genvarname([allvars(v).name,'_clas'])),2));
                dum(1:size(eval(genvarname([allvars(v).name,'_avg'])),1),1:end) = ...
                    eval(genvarname([allvars(v).name,'_avg']));
                eval([genvarname([allvars(v).name,'_clas']),' = cat(3,',...
                    genvarname([allvars(v).name,'_clas']),',' ...
                    'dum);']);
            else
                eval([genvarname([allvars(v).name,'_clas']),' = cat(3,',...
                    genvarname([allvars(v).name,'_clas']),',' ...
                    genvarname([allvars(v).name,'_avg']), ');']);
            end
        end
    end
end

%% k-means clustering
    
% set rates & PAC that could not be detected to 0 
% (to provide a valid range for normalization)
f_Ne_clas(isnan(f_Ne_clas)) = 0;
f_Ni1_clas(isnan(f_Ni1_clas)) = 0;
f_Ni2_clas(isnan(f_Ni2_clas)) = 0;
PAC_clas(isnan(PAC_clas)) = 0;
pBursts_Ne_clas(isnan(pBursts_Ne_clas)) = 0;
pBursts_Ni1_clas(isnan(pBursts_Ni1_clas)) = 0;
pBursts_Ni2_clas(isnan(pBursts_Ni2_clas)) = 0;

vars = {'f_Ne','f_Ni1', 'f_Ni2',...
    'power_LFP_hf','power_LFP_lf',...
    'freq_LFP_lf', 'freq_LFP_hf',...
    'PPC_Ne','PPC_Ni1','PPC_Ni2',...
    'PAC',...
    'pBursts_Ne','pBursts_Ni1','pBursts_Ni2'};

% normalize the variables across all motifs
% convert data to N x p matrix, with N datapoints and p variables
X = [];
for v = 1:length(vars)
    if strcmp(vars{v}(1:3), 'pow')
        dum = eval(['log(',genvarname([vars{v},'_clas']),'(:))']);
    else
        dum = eval([genvarname([vars{v},'_clas']),'(:)']);
    end
    
    % check for division by 0 in dataset
    dum(isinf(dum)) = NaN;
    % standardize
    X(:,v) = (dum-nanmean(dum,1)) ./ nanstd(dum);
end

% set all values that were NaN to 0
delids = find(sum(isnan(X),2)==length(vars));
keepids = setdiff(1:size(X,1),delids);
Xkeep = X(keepids,:);
Xkeep(isnan(Xkeep)) = 0;

% find number of clusters

% k-means clustering
idx = [];
for k = 1:20
    idx(:,k) = kmeans(Xkeep,k,'Distance','sqEuclidean','MaxIter',1000,'Replicates',10);
end

% find optimal number of clusters
va = evalclusters(Xkeep,idx,'CalinskiHarabasz');
figure; plot(va.CriterionValues, 'k'); 
%xlim([4,20])
xlabel('# clusters')
ylabel('Calinski-Harabasz index')

[~,dum] = max(va.CriterionValues(4:end));
nclus = dum+3;

%% Figure 3: plot clusters per motif

sp = [1,2,7:24]; % subplots to use
order = [1,2,...
    3,6,9,12,15,18,...
    4,7,10,13,16,19,...
    5,8,11,14,17,20]; % which motifs should be plotted where
ax = 0:250:5000;

% for 6 clusters (all motifs)
cols = brewermap(8,'Dark2');
dum = brewermap(5,'Paired');
cols = [dum(1,:);cols([3,5,2,6],:)];
dum = brewermap(10,'Blues');
cols = [dum(9,:);cols];
cols = [[1,1,1];cols];

% for 4 clusters (2 motifs)
% cols = brewermap(8,'Dark2');
% cols = [cols([3,5,2],:)];
% dum = brewermap(10,'Blues');
% cols = [cols(1,:);dum(9,:);cols([2,3],:)];
% cols = [[1,1,1];cols];

idx_full = zeros(size(X,1),size(idx,2));
idx_full(keepids,:) = idx;
idx_or = reshape(idx_full(:,nclus), size(eval(genvarname([vars{1},'_clas']))));
figure;
for pp = 1:length(parincl)
    subplot(4,6,sp(pp))
    imagesc(ax,ax,idx_or(:,:,find(parincl == order(pp)))')
    colormap(cols); 
    caxis([0,nclus])
    axis xy
    title(parfiles(order(pp)).name(3:5), 'Interpreter', 'none')
    axis equal; axis tight
    set(gca,'ytick',[0,2500,5000])
    set(gca,'xtick',[0,2500,5000])
end

%% Figure 5: violin plots - plot variable distributions per cluster

vars = {'CV_Ne', 'CV_Ni1','CV_Ni2',...
    'f_Ne','f_Ni1', 'f_Ni2',...
    'f_Ne_norm','f_Ni1_norm', 'f_Ni2_norm',...
    'freq_LFP', 'power_LFP',...
    'power_LFP_hf','power_LFP_lf',...
    'freq_LFP_hf','freq_LFP_lf',... %
    'PPC_Ne','PPC_Ni1','PPC_Ni2',...
    'PAC',...
    'POF_Ne','POF_Ni1','POF_Ni2',...
    'pBursts_Ne','pBursts_Ni1','pBursts_Ni2',...
    };

% convert data to N x p matrix, with N datapoints and p variables
X = [];
for v = 1:length(vars)
    if strcmp(vars{v}(1:3), 'pow')
        dum = eval(['log(',genvarname([vars{v},'_clas']),'(:))']);
    else
        dum = eval([genvarname([vars{v},'_clas']),'(:)']);
    end
    dum(isinf(dum)) = NaN;
    X(:,v) = dum;
end
X(X==0) = NaN;

figure; hold on
varoi = 'freq_LFP';
varnum = find(strcmp(vars, varoi));
vartemp = reshape(X(:,varnum), size(eval(genvarname([vars{1},'_clas']))));
minc = [];
maxc = [];
yax = [0,90];
clusorder = [1,3,4,5,2,6];
for c = 1:nclus
    dum = idx_or == clusorder(c);
    if sum(~isnan(vartemp(dum==1)))< 10
        continue
    end
    scatter(c*ones(sum(dum(:)),1)-0.3+0.2*rand(sum(dum(:)),1), vartemp(dum==1),5,cols(clusorder(c)+1,:))
    violin({vartemp(dum==1)},'x',[c,c+1],'mc','', 'medc','', 'facecolor', cols(clusorder(c)+1,:), 'edgecolor','') 
    scatter(c,nanmean(vartemp(dum==1)),50,cols(clusorder(c)+1,:), 'filled')
minc(c) = min(vartemp(dum==1));
maxc(c) = max(vartemp(dum==1));
end
xlim([0,nclus+1])
xlabel('Cluster ID')
ylabel(varoi, 'Interpreter','none')
ylim(yax);

%% Figure S2: histogram of frequencies 
% run previous section first with varoi = 'freq_LFP';

nbins = 91;
hbars = zeros(nclus,nbins);
bins = linspace(yax(1)+0.5,yax(2)+0.5,nbins);
for c = 1:nclus
    dum = idx_or == clusorder(c);
    if sum(~isnan(vartemp(dum==1)))< 10
        continue
    end
    hbars(c,:) = hist(vartemp(dum==1),bins);
end

figure; hold on
ymax = 600;
freqbound = [4,8,12,30,90];
freqcent = [5.4,9,20,59];
for f = 1:length(freqbound)
   plot([freqbound(f),freqbound(f)],[0,ymax],':k')
end
text(freqcent,0.9*ymax*ones(length(freqcent),1),{'\theta','\alpha','\beta','\gamma'})
b = bar(bins,hbars','stacked','FaceColor','flat', 'EdgeColor','none');
for c = 1:nclus
    b(c).CData = cols(clusorder(c)+1,:);
end
xlim(yax)
ylim([0,ymax])
ylabel('Count')
xlabel(varoi, 'Interpreter','none')

%% Figures 4, 6-8: plot variables - images per motif

vars = {'CV_Ne', 'CV_Ni1','CV_Ni2',...
    'f_Ne','f_Ni1', 'f_Ni2',...
    'f_Ne_norm','f_Ni1_norm', 'f_Ni2_norm',...
    'freq_LFP', 'power_LFP',...
    'power_LFP_hf','power_LFP_lf',...
    'freq_LFP_hf','freq_LFP_lf',... 
    'PPC_Ne','PPC_Ni1','PPC_Ni2',...
    'PAC',...
    'POF_Ne','POF_Ni1','POF_Ni2',...
    'pBursts_Ne','pBursts_Ni1','pBursts_Ni2',...
    };

    PPC_Ne_clas((PPC_Ne_clas)==0) = NaN;
    PPC_Ni1_clas((PPC_Ni1_clas)==0) = NaN;
    PPC_Ni2_clas((PPC_Ni2_clas)==0) = NaN;

% convert data to N x p matrix, with N datapoints and p variables
X = [];
for v = 1:length(vars)
    if strcmp(vars{v}(1:3), 'pow')
        dum = eval(['log(',genvarname([vars{v},'_clas']),'(:))']);
    else
        dum = eval([genvarname([vars{v},'_clas']),'(:)']);
    end
    dum(isinf(dum)) = NaN;
    X(:,v) = dum;
end
delids = find(sum(isnan(X),2)==length(vars));
keepids = setdiff(1:size(X,1),delids);
Xkeep = X(keepids,:);
Xkeep(isnan(Xkeep)) = 0;

dum = brewermap(50,'Spectral');
circcol = [[1,1,1];dum];
normcol = parula(100);
normcol = [[1,1,1];normcol];
figure;
varoi = 'PPC_Ni1';
varnum = find(strcmp(vars, varoi));
vartemp = reshape(X(:,varnum), size(eval(genvarname([vars{1},'_clas']))));
cax = [-0.1,1];
for pp = 1:length(parincl)
    subplot(4,6,sp(pp))
    imagesc(ax,ax,vartemp(:,:,find(parincl == order(pp)))')
    colormap(normcol);  
%     colormap(circcol); % for phases
    axis xy
    title(parfiles(order(pp)).name(3:5), 'Interpreter', 'none')
    caxis(cax)
    set(gca,'ytick',[0,2500,5000])
    set(gca,'xtick',[0,2500,5000])
    axis equal; axis tight
end
subplot(4,6,3); colorbar; caxis(cax)
suptitle(varoi)
