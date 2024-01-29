function [freq, power, spectrum, freqs] = spectralpeak(sig, dt, flim)

% Determines the peak power (power) and corresponding frequency (freq) of a 
% multi-taper spectrogram for signal sig. The spectrum is also returned. 

% Inputs:
% - sig: signal (1D)
% - dt: the time step used for sig (in ms)
% - flim: frequency boundaries for finding the peak

% demean
sig = sig-nanmean(sig);

Fs = (1000/dt); % s-1
% [pxx, freqs] = pmtm(sig,5,length(sig),Fs);
% spectrum = pxx;

SIG = fft(sig);
freqs = Fs/length(sig)*(0:(length(sig)/2));
spectrum = abs(SIG(1:(length(sig)/2+1))/length(sig));


if isempty(flim)
    [~,f1] = min(abs(freqs-1));
    [power, index] = findpeaks(spectrum(f1:end), 'sortstr', 'descend', 'npeaks', 1);
else
    [~,f1] = min(abs(freqs-flim(1)));
    [~,f2] = min(abs(freqs-flim(2)));
    [power, index] = findpeaks(spectrum(f1:f2), 'sortstr', 'descend', 'npeaks', 1);
end

if isempty(index)
    freq = NaN; 
    power = NaN;
else
    freq = freqs(index+f1);
end
