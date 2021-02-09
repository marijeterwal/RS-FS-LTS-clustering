
%{

--- RunSimulations ---
Runs the network simulation defined in the script GetSettings and runs the 
spike- and LFP-based analyses specified in the script SpikeLFPAnalyses.

External dependencies:
- CircStats toolbox

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

clear all

saveID      = 'M1';
saveData    = 1;
saveFigures = 1;
datapath    = 'YourPathHere';

% get parameter settings
% - change GetSettings or replace the line below with any of the specific 
% Motif setting script to run a specific network configuration
GetSettings

% identify loops
duml1 = eval(loop1);    nl1 = length(duml1);
duml2 = eval(loop2);    nl2 = length(duml2);

MemoryAllocation % generate matrices for analysis and saving

for l1 = 1:nl1
    for l2 = 1:nl2
        fprintf('\nRun %i of %i\n', (l1-1)*nl2 + l2, nl1*nl2);

        ChangeParameters % adjust looped variables
        GetConnections % set up connection matrix

        % simulate        
        SimulateNetwork
               
        SpikeLFPAnalyses % compute firing rates, frequencies, PPC, etc.
    end
end

SaveAnalysis