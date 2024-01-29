%{

--- SimulateNetwork ---
This script integrates the spiking neuron network set in RunSimulations and
GetSettings.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

%% initialization

rng('default');
rng(seednr+6);
v = -75*ones(sum(Ntot),1)+10*(rand(sum(Ntot),1)-0.5); % Initial values of v
u = b.*v + d; % Initial values of u
spikes = []; % spike timings
vstore = zeros(sum(Ntot),Nt);
vstore(:,1) = v;

%% inputs

for nr = 1:nregions

    % static input
    Iet{nr} = Ie{nr}*ones(Ne(nr),1) + Ie_std{nr}*randn(Ne(nr),1);
    Ii1t{nr} = Ii1{nr}*ones(Ni1(nr),1) + Ii1_std{nr}*randn(Ni1(nr),1);
    Ii2t{nr} = Ii2{nr}*ones(Ni2(nr),1) + Ii2_std{nr}*randn(Ni2(nr),1);

    % poisson input
    PoissonE = zeros(Ne(nr),Nt/dt); PoissonE(rand(Ne(nr),Nt/dt)<(PRateE{nr}/1000*dt)) = 1;
    PoissonI1 = zeros(Ni1(nr),Nt/dt); PoissonI1(rand(Ni1(nr),Nt/dt)<(PRateI1{nr}/1000*dt)) = 1;
    PoissonI2 = zeros(Ni2(nr),Nt/dt); PoissonI2(rand(Ni2(nr),Nt/dt)<(PRateI2{nr}/1000*dt)) = 1;
    Poisson{nr} = cat(1,PoissonE,PoissonI1,PoissonI2);

    Ip{nr} = zeros(Ntot(nr),1);

    % special thatamic inputs - the same for all neurons within a cell type and region
    % determine which neurons receive what input
    inputid = per{nr}+ramp{nr}+step{nr}+pulse{nr}+multisin{nr};
    if sum(inputid)
        inputselE{nr} = inputid(1)*single(rand(Ne(nr),1)<inputshare{nr});
        inputselI1{nr} = inputid(2)*single(rand(Ni1(nr),1)<inputshare{nr});
        inputselI2{nr} = inputid(3)*single(rand(Ni2(nr),1)<inputshare{nr});
    else
        inputselE{nr} = zeros(Ne(nr),1);
        inputselI1{nr} = zeros(Ni1(nr),1);
        inputselI2{nr} = zeros(Ni2(nr),1);
    end

    rng('default');
    rng(seednr+7); % set seed

    % generate inputs 
    Iper    = per{nr}'*inputamp{nr}*sin(tt*2*pi/Tper{nr});
    Ipulse  = zeros(size(tt)); Ipulse(pulse{nr}==1, inputtime{nr}/dt:(inputtime{nr}+pulselength{nr})/dt) = inputamp{nr};
    Iramp   = ramp{nr}' * 0.01 * tt;
    Istep   = zeros(size(tt)); Istep(step{nr}==1, inputtime{nr}/dt:Nt/dt) = inputamp{nr};
    if ~isempty(multisinloadfile) && sum(multisin{nr})
        load(multisinloadfile)
        Imultisin = multisin{nr}' * (sum(inputamp{nr})+inputamp{nr} * repmat(multisine',[1,length(tt)/length(multisine)]));
    else
        Imultisin = multisin{nr}' * (sum(inputamp{nr})+1/length(sinfreq{nr})*getMultisine(inputamp{nr},sinfreq{nr},sinphase{nr},sinphaseoffset{nr},tt/1000));
    end

    Ithalamus{nr}  = (Iper+Iramp+Istep+Ipulse+Imultisin);
    Pthalamus{nr}   = zeros(3,Nt/dt);
    Pthalamus{nr}(rand(3,Nt/dt)<(Ithalamus{nr}/1000*dt)) = 1;
    It{nr} = zeros(3,1);
end

% r = 0;
% mt = (mu/tau)*dt;
% s2 = sqrt(2*sg*sg*dt/tau);

%% dynamics of synapses

tscaleSyn = [exp(-dt/tauSyn_AMPA)*ones(sum(Ne),1); exp(-dt/tauSyn_GABA_FS)*ones(sum(Ni1),1); exp(-dt/tauSyn_GABA_LTS)*ones(sum(Ni2),1)];
Isyn      = zeros(sum(Ntot),1);
Isyn_delay = zeros(sum(Ntot), DelaySyn/dt);
Isynstore = zeros(nregions,Nt/dt);

%% Run!

rng('default');
rng(seednr+8); % set seed

for t = 2:round(Nt/dt)
    
    for nr = 1:nregions
        % static inputs - noise
        In{nr}    = [Iet{nr} + Ie_std{nr}*randn(Ne(nr),1);...  % noise
                Ii1t{nr} + Ii1_std{nr}*randn(Ni1(nr),1);...
                Ii2t{nr} + Ii2_std{nr}*randn(Ni2(nr),1)];

        % thalamic input
        It{nr} = It{nr}*exp(-dt/tauSyn_AMPA); % mimicking decay of synapses;
        It{nr}(Pthalamus{nr}(:,t)==1) = It{nr}(Pthalamus{nr}(:,t)==1) + 1;
        Itcell{nr} = [inputselE{nr}*It{nr}(1);inputselI1{nr}*It{nr}(2);inputselI2{nr}*It{nr}(3)];

        % poisson input
        Ip{nr} = Ip{nr}*exp(-dt/tauSyn_AMPA); % mimicking decay of synapses
        Ip{nr}(Poisson{nr}(:,t)==1) = Ip{nr}(Poisson{nr}(:,t)==1) + 1;
    end
    % format inputs
    Inall = cat(1,In{:});
    Itall = cat(1,Itcell{:});
    Ipall = cat(1,Ip{:});
        
    % synaptic inputs
    Isyn = Isyn.*tscaleSyn; %with decay - without decay: zeros(Ntot,1);%
    
    % spikes and synapses
    IDfired    = find(v>=0); % indices of spiking neurons
    spikes     = [spikes; t*dt+0*IDfired, IDfired]; 
    v(IDfired) = c(IDfired); % reset
    u(IDfired) = u(IDfired) + d(IDfired); % reset
    Isyn(IDfired) = Isyn(IDfired) + 1;
    
    % store voltage
    vstore(:,t) = v;
    vstore(IDfired,t-1) = 30; % make nice spikes
    
    I = Ipall + Inall + Itall + Syn*Isyn_delay(:,1);% save new input
     
    % store Isyn for delay
    Isyn_delay(:,1:end-1) = Isyn_delay(:,2:end);
    Isyn_delay(:,end) = Isyn;

    for nr = 1:nregions
        Isynstore(nr,t) = sum(abs(Isyn(nctNtot(nr)+1:nctNtot(nr+1))),1);
    end

    v   = v + dt*(0.04*v.^2+5*v+140-u+I);%+dum);
    u   = u + dt*a.*(b.*v-u);  
end

LFP = Isynstore;%mean(vstore,1);

if saveData
    if ~exist([datapath, 'Data/',saveID,'/'])
        mkdir([datapath, 'Data/',saveID,'/'])
    end
    save([datapath, 'Data/',saveID,'/' saveID,'_l1=',num2str(l1),'_l2=',num2str(l2), '_spikeData','.mat'], 'spikes');
    save([datapath, 'Data/',saveID,'/' saveID,'_l1=',num2str(l1),'_l2=',num2str(l2), '_LFP','.mat'], 'LFP','Isynstore','Ithalamus');
end
    