% Motif XX

%% loops
loop1 = 'inputamp{1}';%'PRateI1'; 
loop2 = 'NaN';%'PRateE';

%% time
Nt      = 1500; %ms
dt      = 0.2;
dthist  = 2; %ms
tsel    = dt;%300; %ms - Time after onset to be excluded from analyses
tt      = 0:dt:Nt-dt;

% # of neurons
Ne      = 800;%[800,800]; % RS       
Ni1     = 200;%[100,200]; % FS     
Ni2     = 00;%[100,0]; % LTS

Ni      = Ni1+Ni2;
Ntot    = Ne+Ni;
nregions = length(Ne);

%% random numbers

seednr = 1; % seed

rng('default')
rng(seednr(1)+1);     re      = rand(sum(Ne),1);          
rng(seednr(1)+2);     ri1     = rand(sum(Ni1),1);    
rng(seednr(1)+3);     ri2     = rand(sum(Ni2),1);

%% neuron parameters

nctNtot = [0 cumsum(Ntot)];
nctNe = [0 cumsum(Ne)];
nctNi1 = [0 cumsum(Ni1)];
nctNi2 = [0 cumsum(Ni2)];

[a,b,c,d] = deal([]);
for nr = 1:length(Ne)
    a       =   [a; [0.02*ones(Ne(nr),1);
                0.1+0.08*ri1(nctNi1(nr)+1:nctNi1(nr+1));
                0.02+0.005*ri2(nctNi2(nr)+1:nctNi2(nr+1))]];
    b       =   [b; [0.2*ones(Ne(nr),1);
                0.2-0.05*ri1(nctNi1(nr)+1:nctNi1(nr+1));
                0.25-0.05*ri2(nctNi2(nr)+1:nctNi2(nr+1))]];
    c       =   [c; [-65+15*re(nctNe(nr)+1:nctNe(nr+1)).^2;
                -65*ones(Ni(nr),1)]];
    d       =   [d; [8-6*re(nctNe(nr)+1:nctNe(nr+1)).^2;
                2*ones(Ni(nr),1)]];
end

%% Inputs
Ie      = {0,0};
Ii1     = {0,0};
Ii2     = {0,0};

Ie_std  = {2,2};
Ii1_std = {2,2};
Ii2_std = {2,2};

PRateE = {1500,500};%0:500:2500;
PRateI1 = {1000,500};%0:500:2500;
PRateI2 = {100,500};%1000;

per    = {[0,0,0],[0,0,0]};
ramp   = {[0,0,0],[0,0,0]};
step   = {[0,0,0],[0,0,0]};
pulse  = {[1,0,0],[0,0,0]};
multisin  = {[0,0,0],[0,0,0]};

inputshare = {1,0};
inputtime = {705,710};
inputamp = {500,0};%1000;%0:100:1000;
Tper    = {100,100};%20;
pulselength = {50,0};
sinfreq = {[1:2:15,19,23],[]};%[3,7,13];%
sinphase = {[],[]};%[0,pi/3,-pi/3];%
sinphaseoffset = {0*pi,0};
multisinloadfile = 'U:\DBI2 BME\User Marije\Multisines\Vlaar_multisine_fs5000_4.mat';

% stimtype = 'fullfield';
% sg = 1;
% mu = 1;
% tau = 20; %ms

%% synapses

AvgWeights = 2;
StdWeights = 1;

% SynConMat = [0.5 -1 0; % weights RS-FS - motif I
%              0.5 -1 0;
%              0  0  0];

% SynConMat = [0.5 -1 -1; % weights RS-FS-LTS - motif VII or VIII
%              0.5 -1 -1;
%              0  0  0];

SynConMat = {};
SynConMat{1,1} =    [0.5 -1 -1; % weights RS-FS-LTS - motif XIX or XX
                    0.5 -1 -1;
                    0.5  0  0];
SynConMat{2,2} =    [0.5 -1 -1; % weights RS-FS-LTS - motif XIX or XX
                    0.5 -1 -1;
                    0.5  0  0];
SynConMat{1,2} =    [0.5 -1 0; % weights
                     0.5 0 0;
                     0.5 0 0];
SynConMat{2,1} =    [0.5 -1 0; % weights
                     0.5 0 0;
                     0.5 0 0];

pconSyn = {};
pconSyn{1,1} =  [0.05 0.30 0.40; % RS-FS-LTS
                0.10 0.30 0.20;
                0.10 0.20 0.00];
pconSyn{2,2} =  [0.05 0.30 0.40; % RS-FS-LTS
                0.10 0.30 0.20;
                0.10 0.20 0.00];
pconSyn{1,2} = [0.05 0.0 0.0; % RS-FS-LTS
                0.0 0.0 0.0;
                0.0 0.0 0.0];
pconSyn{2,1} = [0.05 0.0 0.0; % RS-FS-LTS
                0.0 0.0 0.0;
                0.0 0.0 0.0];

% time scale of synaptic decay
tauSyn_GABA_FS = 3; %ms
tauSyn_GABA_LTS = 6; %ms
tauSyn_AMPA = 2; %ms

% synaptic delay
DelaySyn = 1; %ms

%% analysis

% gaussian smoothing
sigma = 50; %ms

burstISI = 10; %ms

%% Save Parameters & Settings

if exist([datapath, '/', 'Figures/', saveID, '/'], 'dir') == 0
    mkdir([datapath, 'Data/', saveID, '/']); 
    mkdir([datapath, 'Figures/', saveID, '/']);    
    mkdir([datapath, 'Analysis/', saveID, '/']);    
end

if saveData || saveFigures
    save([datapath, 'Settings/', saveID, '_Settings','.mat']);
    Svars = who;
    fileID = fopen([datapath, 'Settings/', saveID, '_Settings','.txt'], 'wt+');
    for s = 1:length(Svars)
        if isequal(Svars{s}, 'Svars')
            continue
        end
        try
            fprintf(fileID, [Svars{s} ' = \r\n' mat2str(eval(Svars{s})) '\r\n\r\n']); %
        end
    end
    fclose(fileID);
end