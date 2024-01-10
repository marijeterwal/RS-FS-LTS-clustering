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
v = -75*ones(Ne+Ni,1)+10*(rand(Ne+Ni,1)-0.5); % Initial values of v
u = b.*v + d; % Initial values of u
spikes = []; % spike timings
vstore = zeros(Ntot,Nt);
vstore(:,1) = v;

%% inputs

Iet = Ie*ones(Ne,1) + Ie_std*randn(Ne,1);
Ii1t = Ii1*ones(Ni1,1) + Ii1_std*randn(Ni1,1);
Ii2t = Ii2*ones(Ni2,1) + Ii2_std*randn(Ni2,1);

PoissonE = zeros(Ne,Nt/dt); PoissonE(rand(Ne,Nt/dt)<(PRateE/1000*dt)) = 1;
PoissonI1 = zeros(Ni1,Nt/dt); PoissonI1(rand(Ni1,Nt/dt)<(PRateI1/1000*dt)) = 1;
PoissonI2 = zeros(Ni2,Nt/dt); PoissonI2(rand(Ni2,Nt/dt)<(PRateI2/1000*dt)) = 1;
Poisson = cat(1,PoissonE,PoissonI1,PoissonI2);

Ip = zeros(Ntot,1);

r = 0;
mt = (mu/tau)*dt;
s2 = sqrt(2*sg*sg*dt/tau);

if sum(inputid)
    inputselE = inputid(1)*single(rand(Ne,1)<inputshare);
    inputselI1 = inputid(2)*single(rand(Ni1,1)<inputshare);
    inputselI2 = inputid(3)*single(rand(Ni2,1)<inputshare);
else
    inputselE = zeros(Ne,1);
    inputselI1 = zeros(Ni1,1);
    inputselI2 = zeros(Ni2,1);
end

rng('default');
rng(seednr+7); % set seed

% special currents
Iper    = per'*inputamp*sin(tt*2*pi/Tper);
Ipulse  = pulse' * zeros(size(tt)); Ipulse(pulse==1, inputtime/dt:(inputtime+pulselength)/dt) = inputamp;
Iramp   = ramp' * 0.01 * tt;
Istep   = step'* zeros(size(tt)); Istep(step==1, inputtime/dt:Nt/dt) = inputamp;
Imultisin = multisin' * (sum(inputamp)+1/length(sinfreq)*getMultisine(inputamp,sinfreq,sinphase,tt/1000));

Ithalamus   = (Iper+Iramp+Istep+Ipulse+Imultisin); 
Pthalamus   = zeros(3,Nt/dt);
Pthalamus(rand(3,Nt/dt)<(Ithalamus/1000*dt)) = 1;
It = zeros(3,1);


%% dynamics of synapses

tscaleSyn = [exp(-dt/tauSyn_AMPA)*ones(Ne,1); exp(-dt/tauSyn_GABA_FS)*ones(Ni1,1); exp(-dt/tauSyn_GABA_LTS)*ones(Ni2,1)];
Isyn      = zeros(Ntot,1);
Isyn_delay = zeros(Ntot, DelaySyn/dt);
Isynstore = zeros(1,Nt/dt);

%% Run!

rng('default');
rng(seednr+8); % set seed

for t = 2:round(Nt/dt)
    
    % static inputs - noise
    In    = [Iet + Ie_std*randn(Ne,1);...  % noise 
            Ii1t + Ii1_std*randn(Ni1,1);...
            Ii2t + Ii2_std*randn(Ni2,1)];
    
    % thalamic input
    It = It*exp(-dt/tauSyn_AMPA); % mimicking decay of synapses;
    It(Pthalamus(:,t)==1) = It(Pthalamus(:,t)==1) + 1;   

    % poisson input    
    Ip = Ip*exp(-dt/tauSyn_AMPA); % mimicking decay of synapses
    Ip(Poisson(:,t)==1) = Ip(Poisson(:,t)==1) + 1;
    
        
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
    
    I = Ip + Syn*Isyn_delay(:,1) + In + [inputselE*It(1);inputselI1*It(2);inputselI2*It(3)];% save new input
     
    % store Isyn for delay
    Isyn_delay(:,1:end-1) = Isyn_delay(:,2:end);
    Isyn_delay(:,end) = Isyn;

    Isynstore(t) = sum(abs(Isyn),1);

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
    