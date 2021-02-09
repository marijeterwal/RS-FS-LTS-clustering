%{

--- GetSettings ---
This script defines all parameters for the network simulations run by RunSimulations
and analyses in SpikeLFPAnalyses.

If they do not exist, it makes directories for storing data, analysis
results, figures and settings. It stores a .mat file with the current workspace
and prints all settings to a .txt file for future reference.

This script is currently set for fI curve simulations (Figure 1) for the RS 
cell population across 10 seeds. Connections have been turned off. 

Requires the following variables to be defined in the workspace:
- 'datapath': string with the path to the main directory used for storing the simulation data
- 'saveID': string used as a unique identified for the simulation

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

%% parameter loop

% select the parameters to loop over
% for no parameter loop set to 'NaN'
loop1 = 'PRateE'; 
loop2 = 'seednr';

%% time

Nt      = 2300; %ms
dt      = 0.2;
tsel    = 300; %ms - Time after onset to be excluded from analyses
tt      = 0:dt:Nt-dt;

% # of neurons
Ne      = 800;% RS       
Ni1     = 100; % FS     
Ni2     = 100; % LTS

Ni      = Ni1+Ni2;
Ntot    = Ne+Ni;

%% Random numbers

seednr = 1:3; % seed

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
            2*ones(Ni1,1);
            2*ones(Ni1,1)];

%% Inputs

% static currents
Ie      = 0;
Ii1     = 0;
Ii2     = 0;

% noise
Ie_std  = 1;
Ii1_std = 1;
Ii2_std = 1;

PRateE = 0:1000:5000;
PRateI1 = 0;
PRateI2 = 0;

% special inputs per cell type:
per    = [0,0,0]; % periodic
ramp   = [0,0,0]; % ramp
step   = [0,0,0]; % step
pulse  = [0,0,0]; % pulse
inputid = per+ramp+step+pulse;

inputshare = 1;
inputtime = 1000;
inputamp = 5;
Tper    = 20;
pulselength = 200;

stimtype = 'fullfield';
sg = 1;
mu = 1;
tau = 20; %ms

%% Synapses

AvgWeights = 2;
StdWeights = 1;

SynConMat = zeros(3);
% SynConMat = [0.5 -1 -1; % weights RS-FS-LTS
%              0.5 -1 -1;
%              0.5  -1  0];

pconSyn = [0.05 0.30 0.40; % RS-FS-LTS
           0.10 0.30 0.20;
           0.10 0.20 0.00];

% time scale of synaptic decay
tauSyn_GABA_FS = 3; %ms
tauSyn_GABA_LTS = 6; %ms
tauSyn_AMPA = 2; %ms

DelaySyn = 1; %ms

%% Analysis

% gaussian smoothing
sigma = 50; %ms

burstISI = 10; %ms

%% Save Parameters & Settings

if exist([datapath, 'Figures/', saveID, '/'], 'dir') == 0
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
        fprintf(fileID, [Svars{s} ' = \r\n' mat2str(eval(Svars{s})) '\r\n\r\n']); %
    end
    fclose(fileID);
end