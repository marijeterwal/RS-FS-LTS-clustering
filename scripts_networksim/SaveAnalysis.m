%{

--- SaveAnalysis ---
Saves all variables defined in MemoryAllocation.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}
if saveData
    save([datapath, 'Analysis/', saveID, '_Analysis.mat'], varNames{:}); % varNames from MemoryAllocation script
end