
%{

--- PlotfIcurve ---
Loads the firing rate data from a specified simulation run and plots the
firing rate distribution against the loop1 variable, which should be input 
current or Poisson input rate (see GetSettings). 
The loop2 variable can be seednr, or NaN.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

clear all
clc

saveID = 'M1'; % run identifier
datapath = 'YourPathHere';
saveFigures = false;

celltype = 'Ne'; % cell type to plot;
varname = strcat('fcell_', celltype);

%% load data

load([datapath,'Analysis\',saveID,'_Analysis.mat'])
load([datapath,'Settings\',saveID,'_Settings.mat'])

fr = eval(['reshape(permute(',varname,',[2,1,3]),[size(',varname,',2),',num2str(length(seednr)),'*',celltype,'])']);

xax = eval(loop1);
xaxlabel = loop1;

yaxlabel = 'Firing rate (Hz)';

%% plot f-I curve

fr_avg = mean(fr,2);
fr_std = std(fr,0,2);
fr_low = quantile(fr,0.025,2);
fr_high = quantile(fr,0.975,2);

figure; hold on
fill([xax,fliplr(xax)],[fr_high;flipud(fr_low)],[0.6,0.6,0.6],...
    'linestyle','none');
plot(xax, fr_avg, 'k', 'linewidth', 2);
xlim([xax(1),xax(end)]);
ylim([-2, 120]);
xlabel(xaxlabel);
ylabel(yaxlabel);

%% save figure
if saveFigures
    saveas(gcf,[datapath, 'Figures/', saveID, '_fIcurve'], 'fig');
end
