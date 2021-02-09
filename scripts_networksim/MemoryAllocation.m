%{

--- MemoryAllocation ---
Memory allocation for all variables in SpikeLFPAnalyses. 
Only variables defined below are saved by SaveAnalysis.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

L1 = length(duml1);
L2 = length(duml2);

% L1 x L2 matrices
varNames1 = {'f_Ne', 'f_Ni1', 'f_Ni2',...
    'freq_LFP', 'power_LFP', 'freq_LFP_hf', 'power_LFP_hf', 'freq_LFP_lf', 'power_LFP_lf',...
    'PAC',...
    'PPC_Ne', 'PPC_Ni1', 'PPC_Ni2','PPC_all',...
    'pBursts_Ne', 'pBursts_Ni1', 'pBursts_Ni2',...
    'POF_Ne', 'POF_Ni1','POF_Ni2',...
    'Rf_Ne', 'Rf_Ni1','Rf_Ni2',...
    'pvalRf_Ne', 'pvalRf_Ni1','pvalRf_Ni2'};

for i = 1:length(varNames1)
    eval([genvarname(varNames1{i}) '= zeros(L1,L2);']);
end

% L1 x L2 cells
varNames2 = {'spectrum_LFP',...
    'phist_Ne','phist_Ni1', 'phist_Ni2'};

for i = 1:length(varNames2)
    eval([genvarname(varNames2{i}) '= cell(L1,L2);']);
end

varNames3 = {'fcell_Ne','fcell_Ni1','fcell_Ni2'};
for i = 1:length(varNames3)
    eval([genvarname(varNames3{i}) '= [];']);
end

varNames = cat(2,varNames1,varNames2,varNames3);

