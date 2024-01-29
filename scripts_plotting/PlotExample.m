%{

--- PlotExample ---
Loads the spike and LFP data from a specified simulation run and plots the
raster plot, LFP trace and spectrum of the LFP.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

% for loading
datapath = 'U:\DBI2 BME\User Marije\Simulations\Data\';
run = 'S19'; 
l1 = '1'; % in case loops were used in this simulation, specify which value to use
l2 = '1';

% for plotting:
tlim = [500,1200];

cl = [];
cl.red = [190,30,45]/255;
cl.blue = [42,56,143]/255;
cl.green = [0,103,56]/255;
cl.purple = [102, 45,145]/255;
cl.orange = [247,148,29]/255;

regcol = [cl.purple; cl.orange];

%% load
if  exist([datapath, 'Data\',run,'\',run,'_l1=',l1,'_l2=',l2,'_LFP.mat'],'file')
    load([datapath, 'Data\',run,'\',run,'_l1=',l1,'_l2=',l2,'_LFP.mat'])
    load([datapath, 'Data\',run,'\',run,'_l1=',l1,'_l2=',l2,'_spikeData.mat'])
    load([datapath,'Settings\',run,'_Settings.mat'])
else
    error('Data not found')
end

%% plot

figure(1);
subplot(nregions+1,1,1); hold on
for nr = 1:nregions
    plot(tt, LFP(nr,:), 'color',regcol(nr,:));
end
ylabel('LFP (a.u.)')
ylim([0,100])% ylim([-130,-40])
xlim(tlim);
legend({'Circuit # 1','Circuit #2'})

for nr = 1:nregions
    subplot(nregions+1,1,nr+1); hold on
    regids = spikes(:,2)>=nctNtot(nr)+1 & spikes(:,2)<nctNtot(nr+1);
    data = spikes(regids==1,:);
    data(:,2) = data(:,2) - nctNtot(nr);
    plot(data(data(:,2) > Ne(nr)+Ni1(nr),1),data(data(:,2) > Ne(nr)+Ni1(nr),2),'.','color',cl.green);
    plot(data(data(:,2) > Ne(nr) & data(:,2) <= Ne(nr)+Ni1(nr),1),data(data(:,2) > Ne(nr) & data(:,2) <= Ne(nr)+Ni1(nr),2),'.','color',cl.blue);
    plot(data(data(:,2) <= Ne(nr),1),data(data(:,2) <= Ne(nr),2),'.','color',cl.red);
    ylabel('Neuron number'); xlabel('Time (ms)')
    plot(tt, Ithalamus{nr}(1,:),'k')
    ylim([0 Ntot(nr)]);
    xlim([tlim]);
end

%% spectra

f2 = figure(2); hold on
f3 = figure(3); hold on

for nr = 1:nregions
    LFP_sel = LFP(nr,tsel/dt:end);
    [~,~,spectrum,freqs] = spectralpeak(LFP_sel, dt,[2,150]);
    [~,~,spectrumT,freqsT] = spectralpeak(Ithalamus{nr}(1,tsel/dt:end)/10, dt,[2,150]);

    figure(f2)
    semilogy(freqsT,spectrumT,'color',[0.6,0.6,0.6])
    semilogy(freqs,spectrum,'color',regcol(nr,:),'linewidth',2)
    xlim([0.5,100])
    ylim([0,5])
    ylabel('Power (a.u.)')
    xlabel('Frequency (kHz)')

    figure(f3);
    subplot(nregions,1,nr)
    [spec,freq] = cwt(LFP(nr,:),1000/dt);
    pcolor(tt,freq,abs(spec)); shading flat
%     set(gca,'YScale','log')
    colorbar
    ylim([3,100])
    ylabel('Frequency (Hz)')
    xlabel('Time (ms)')
end

%% plot input

figure;
plot(tt, cat(1,Ithalamus{:}),'k')
ylabel('Thalamic input (Hz)')
ylim([0,1000])
xlim(tlim);