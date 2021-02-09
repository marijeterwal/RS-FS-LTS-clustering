%{

--- PlotExample ---
Loads the spike and LFP data from a specified simulation run and plots the
raster plot, LFP trace and spectrum of the LFP.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

% for loading
datapath = 'YourPathHere';
run = 'M1'; 
l1 = '1'; % in case loops were used in this simulation, specify which value to use
l2 = '1';

% for plotting:
tlim = [500,1500];

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
ylim([-130,-40])
xlim(tlim);

subplot(212); hold on
plot(spikes(spikes(:,2) > Ne+Ni1,1),spikes(spikes(:,2) > Ne+Ni1,2),'.','color',cl.green);
plot(spikes(spikes(:,2) > Ne & spikes(:,2) <= Ne+Ni1,1),spikes(spikes(:,2) > Ne & spikes(:,2) <= Ne+Ni1,2),'.','color',cl.blue);
plot(spikes(spikes(:,2) <= Ne,1),spikes(spikes(:,2) <= Ne,2),'.','color',cl.red);
ylabel('Neuron number'); xlabel('Time (ms)')
ylim([0, Ntot+1]);
xlim([tlim]);

%% spectra

 LFP_sel = LFP(tsel/dt:end);
[~,~,spectrum,freqs] = spectralpeak(LFP_sel, dt,[2,150]);

figure;
plot(freqs,spectrum,'k')
xlim([2,100])
ylim([0,5])
ylabel('Power (a.u.)')
xlabel('Frequency (Hz)')


