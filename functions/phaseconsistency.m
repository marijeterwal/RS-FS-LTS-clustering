
function [ppc] = phaseconsistency(spikes, ref, dt)

% Determines the Pairwise Phase Consistency (PPC)
% see Vinck et al., Journal of Computational Neuroscience, 2012

% Inputs:
% - spikes: variable of size nspikes x 2, with the first column containing
% the spike times (in ms), and the second column containing the neuron ID.
% - ref: phase trace (in radians)
% - dt: the time step used for ref (in ms)

if size(spikes,1) > 150000
    warning('Spikes is large, PPC will take a lot of time')
end

spikeTimes = round(spikes(:,1)/dt); % time IDs of spikes
if isempty(spikeTimes)
    warning('No spikes found - no PPC computed')
end
phaseNr = ref(spikeTimes); % phases
Ns = length(spikeTimes); % number of spikes

U1 = cos(phaseNr);
U2 = sin(phaseNr);
U = [U1',U2'];

Pl_n = zeros(Ns,1);
for ns = 1:Ns % spikes
    Pl_n = Pl_n + U * [U1(ns),U2(ns)]';
end
ppc = (sum(Pl_n(:))-Ns)/(Ns*(Ns-1));

end