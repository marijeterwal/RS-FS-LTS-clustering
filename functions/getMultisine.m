function sigSin = getMultisine(amps,freqs,phases, phaseoffset, times)

if ~isempty(phases) && length(freqs) ~= length(phases)
    error('Number of frequencies does not match number of phases')
elseif isempty(phases)
    warning('Random phases selected')
    phases = 2*pi*rand(length(freqs),1)-pi;
end

if length(amps)==1 && length(amps) ~= length(freqs)
    warning('Only one amplitude defined, using the same amp for all frequencies')
    amps = amps*ones(size(freqs));
elseif length(amps) ~= length(freqs)
    error('Number of amplitudes does not match number of frequencies')
end

sig = zeros(length(amps),length(times));
for a = 1:length(amps)
    sig(a,:) = amps(a)*sin(times*2*pi*freqs(a)-phases(a)+phaseoffset);
end
sigSin = sum(sig,1);

end