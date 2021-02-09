% Motif XIX

%% loops
loop1 = 'PRateI1'; 
loop2 = 'PRateE';

%% time
Nt      = 2300; %ms
dt      = 0.2;
dthist  = 2; %ms
tsel    = 300; %ms - Time after onset to be excluded from analyses
tt      = 0:dt:Nt-dt;

% # of neurons
Ne      = 800; % RS       
Ni1     = 100; % FS     
Ni2     = 100; % LTS

Ni      = Ni1+Ni2;
Ntot    = Ne+Ni;

%% Random numbers

seednr = 1; % seed

rng('default')
rng(seednr(1)+1);     re      = rand(Ne,1);          
rng(seednr(1)+2);     ri1     = rand(Ni1,1);    
rng(seednr(1)+3);     ri2     = rand(Ni2,1);

%% neuron parameters
a       =   [0.02*ones(Ne,1);
            0.1+0.08*ri1;
            0.02+0.005*ri2];
b       =   [0.2*ones(Ne,1);
            0.2-0.05*ri1;
            0.25-0.05*ri2];
c       =   [-65+15*re.^2;
            -65*ones(Ni,1)];
d       =   [8-6*re.^2;
            2*ones(Ni,1)];

%% Inputs
Ie      = 0;
Ii1     = 0;
Ii2     = 0; 

Ie_std  = 1;
Ii1_std = 1;
Ii2_std = 1;

PRateE = 0:250:5000;
PRateI1 = 0:250:5000;
PRateI2 = 0;

per    = [0,0,0];
ramp   = [0,0,0];
step   = [0,0,0];
pulse  = [0,0,0];
inputid = per+ramp+step+pulse;

inputshare = 0;
inputtime = 100;
inputamp = 10;
Tper    = 20; 

stimtype = 'fullfield';
sg = 1;
mu = 1;
tau = 20; %ms

%% Synapses

AvgWeights = 2;
StdWeights = 1;

SynConMat = [0.5 -1 -1; % weights RS-FS-LTS
             0.5 -1 -1;
             0.5  0  0];

pconSyn = [0.05 0.30 0.40; % RS-FS-LTS
           0.10 0.30 0.20;
           0.10 0.20 0.00];

% time scale of synaptic decay
tauSyn_GABA_FS = 3; %ms
tauSyn_GABA_LTS = 6; %ms
tauSyn_AMPA = 2; %ms

% synaptic delay
DelaySyn = 1; %ms

%% Analysis

% gaussian smoothing
sigma = 50; %ms

burstISI = 10; %ms

%% Save Parameters & Settings

if exist([datapath, '/', 'Figures/', saveID, '/'], 'dir') == 0
    mkdir([datapath, 'Data/', saveID, '/']); 
    mkdir([datapath, 'Figures/', saveID, '/']);    
    mkdir([datapath, 'Settings/', saveID, '/']);  
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