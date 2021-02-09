
function [pof,phist,r,pval] = phaseoffiring(spikes, ref, dt, bins)

% Determines the mean phase of firing (pof) relative to a reference trace 
% and outputs a phase histogram of all spikes (phist), Rayleigh's Z (r) and 
% corresponding p-value (pval). 

% Inputs:
% - spikes: variable of size nspikes x 2, with the first column containing
% the spike times (in ms), and the second column containing the neuron ID.
% - ref: phase trace (in radians)
% - dt: the time step used for ref (in ms)
% - bins: number of bins or bin edges to use for histogram

spikeTimes = round(spikes(:,1)/dt); % time IDs of spikes
if isempty(spikeTimes)
    warning('No spikes found - no phase analysis')
    pof = NaN;
    phist = [];
    r = NaN;
    pval = NaN;
    return
end
phaseNr = ref(spikeTimes); % phases

% mean phase of firing
pof = angle(sum(exp(1i*phaseNr)));

% firing histogram
phist = hist(phaseNr,bins);

try
    % rayleigh's test
    [pval, r] = circ_rtest(phaseNr);
catch
    warning('Circular Statistics toolbox not found, skipping...')
    pval = NaN; 
    r = NaN;
end

end