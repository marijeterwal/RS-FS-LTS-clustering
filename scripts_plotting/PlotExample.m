%{

--- PlotExample ---
Loads the spike and LFP data from a specified simulation run and plots the
raster plot, LFP trace and spectrum of the LFP.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

% for loading
datapath = 'U:\4DEEG\User Marije\Simulations\Data\';
run = 'M14'; 
l1 = '1'; % in case loops were used in this simulation, specify which value to use
l2 = '1';

% for plotting:
tlim = [000,2300];

cl = [];
cl.red = [190,30,45]/255;
cl.blue = [42,56,143]/255;
cl.green = [0,103,56]/255;

%% load
if  exist([datapath, 'Data\',run,'\',run,'_l1=',l1,'_l2=',l2,'_LFP.mat'],'file')
    load([datapath, 'Data\',run,'\',run,'_l1=',l1,'_l2=',l2,'_LFP.mat'])
    load([datapath, 'Data\',run,'\',run,'_l1=',l1,'_l2=',l2,'_spikeData.mat'])
    load([datapath,'Settings\',saveID,'_Settings.mat'])
else
    error('Data not found')
end

%% plot

figure;
subplot(211); hold on
plot(tt, LFP, 'k');
ylabel('LFP (a.u.)')
ylim([0,100])% ylim([-130,-40])
xlim(tlim);

subplot(212); hold on
plot(spikes(spikes(:,2) > Ne+Ni1,1),spikes(spikes(:,2) > Ne+Ni1,2),'.','color',cl.green);
plot(spikes(spikes(:,2) > Ne & spikes(:,2) <= Ne+Ni1,1),spikes(spikes(:,2) > Ne & spikes(:,2) <= Ne+Ni1,2),'.','color',cl.blue);
plot(spikes(spikes(:,2) <= Ne,1),spikes(spikes(:,2) <= Ne,2),'.','color',cl.red);
ylabel('Neuron number'); xlabel('Time (ms)')
plot(tt, Ithalamus(1,:),'k')
ylim([0, Ntot+1]);
xlim([tlim]);

%% spectra

LFP_sel = LFP(tsel/dt:end);
[~,~,spectrum,freqs] = spectralpeak(LFP_sel, dt,[2,150]);
[~,~,spectrumT,freqsT] = spectralpeak(Ithalamus(1,tsel/dt:end)/10, dt,[2,150]);

figure;hold on
semilogy(freqs,spectrum,'k')
semilogy(freqsT,spectrumT,'r')
xlim([0.5,100])
ylim([0,5])
ylabel('Power (a.u.)')
xlabel('Frequency (kHz)')

figure; 
cwt(LFP,1000/dt)
ylabel('Frequency (Hz)')
xlabel('Time (ms)')

%% plot input

figure;
plot(tt, Ithalamus(1,:),'k')
ylabel('Thalamic input (Hz)')
ylim([0,1000])
xlim(tlim);