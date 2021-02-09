%{

--- SpikeLFPAnalyses ---
This script deals with all the spike- and LFP-based analyses needed for the
clustering algorithm in ClusterAnalysis. Memory for the variables filled
here is allocated in MemoryAllocation and variables are saved by
SaveAnalysis.

External dependencies:
- CircStats toolbox

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

if ~isempty(spikes)
    
    data = spikes(spikes(:,1) > tsel,:);
    
    data_Ne = data(data(:,2) <= Ne,:);
    data_Ni1 = data(data(:,2) > Ne & data(:,2)<= Ne+Ni1,:);
    data_Ni2 = data(data(:,2) > Ne+Ni1,:);
    
    %% firing rate
    
    f_Ne(l1,l2) = (length(data_Ne)/Ne) * (1000/(Nt-tsel));
    f_Ni1(l1,l2) = (length(data_Ni1)/Ni1) * (1000/(Nt-tsel));
    f_Ni2(l1,l2) = (size(data_Ni2,1)/Ni2) * (1000/(Nt-tsel));
    
    % firing rate per cell
    for n = 1:Ne
        fcell_Ne(n,l1,l2) = length(data_Ne(data_Ne(:,2)==n,1)) * (1000/(Nt-tsel));
    end
    for n = 1:Ni1
        fcell_Ni1(n,l1,l2) = length(data_Ni1(data_Ni1(:,2)==n,1)) * (1000/(Nt-tsel));
    end
    for n = 1:Ni2
        fcell_Ni2(n,l1,l2) = length(data_Ni2(data_Ni2(:,2)==n,1)) * (1000/(Nt-tsel));
    end
    
    %% LFP freq and phase
    
    LFP_sel = LFP(tsel/dt+1:end);
    [freq_LFP(l1,l2), power_LFP(l1,l2), spectrum_LFP{l1,l2},freqs] = spectralpeak(LFP_sel, dt,[2,150]);
    
    % low freq peak
    [freq_LFP_lf(l1,l2), power_LFP_lf(l1,l2), ~] = spectralpeak(LFP_sel, dt,[2,30]);
    % gamma peak
    [freq_LFP_hf(l1,l2), power_LFP_hf(l1,l2), ~] = spectralpeak(LFP_sel, dt,[30,150]);
    
    %% phase amplitude coupling
     
    % exclude 1st harmonic
	dum = freq_LFP_hf(l1,l2) / freq_LFP_lf(l1,l2);
    if abs(dum-2) < 0.1
        freq_LFP_hf(l1,l2) = NaN;
        power_LFP_hf(l1,l2) = NaN;
        PAC(l1,l2) = NaN;
    elseif freq_LFP_hf(l1,l2) < 40 || power_LFP_lf(l1,l2)<1 || power_LFP_hf(l1,l2)<1
        PAC(l1,l2)= NaN;
    else
        % phase
        [bfilt,afilt] = butter(2, [2,30]/(500/dt));
        filtLFP = filtfilt(bfilt,afilt,LFP_sel);
        hlbrt_lf = hilbert(filtLFP-mean(filtLFP));
        hlbrt_lf = hlbrt_lf ./ norm(hlbrt_lf);

        % amplitude
        [bfilt,afilt] = butter(4, [30,150]/(500/dt));
        filtLFP = filtfilt(bfilt,afilt,LFP_sel);
        hlbrt_hf = abs(hilbert(filtLFP-mean(filtLFP)));
        hlbrt_hf = hlbrt_hf ./ norm(hlbrt_hf);

        PAC(l1,l2) = abs(hlbrt_hf * hlbrt_lf');
    end
    
    %% instantaneous phase

    dumfr = freq_LFP(l1,l2);

    if ~isnan(dumfr)
        
        [bfilt,afilt] = butter(2, [max(2,dumfr-5),dumfr+5]/(500/dt));
        filtLFP = filtfilt(bfilt,afilt,LFP);
        hlbrt = hilbert(filtLFP-mean(filtLFP));
        instphase = angle(hlbrt);
        
        tgs = -100:dt:100;
        gs = dt/(sigma*sqrt(2*pi))*exp(-tgs.^2/(2*sigma^2));
        instfreq = conv((1000/dt)/(2*pi)*diff(unwrap(angle(hlbrt))),gs, 'same');
        
        %% pairwise phase consistency
        
        reference = instphase;
        if ~isempty(data_Ne)
            dumdata = data_Ne(data_Ne(:,2) <= 100,:);
            [PPC_Ne(l1,l2)] = phaseconsistency(dumdata, reference, dt);
        end
        if ~isempty(data_Ni1)
            [PPC_Ni1(l1,l2)] = phaseconsistency(data_Ni1, reference, dt);
        end
        if ~isempty(data_Ni2)
            [PPC_Ni2(l1,l2)] = phaseconsistency(data_Ni2, reference, dt);
        end
        if ~isempty(data)
            [PPC_all(l1,l2)] = phaseconsistency(data, reference, dt);
        end
        
        %% phase of firing
        bins = linspace(-pi,pi,15);
        [POF_Ne(l1,l2) phist_Ne{l1,l2},Rf_Ne(l1,l2),pvalRf_Ne(l1,l2)] = phaseoffiring(data_Ne, reference,dt, bins);
        [POF_Ni1(l1,l2) phist_Ni1{l1,l2},Rf_Ni1(l1,l2),pvalRf_Ni1(l1,l2)] = phaseoffiring(data_Ni1, reference,dt, bins);
        [POF_Ni2(l1,l2) phist_Ni2{l1,l2},Rf_Ni2(l1,l2),pvalRf_Ni2(l1,l2)] = phaseoffiring(data_Ni2, reference,dt, bins);
        
    end
    
    %% bursts
    nsingle = zeros(Ntot,1);
    nburst = zeros(Ntot,1);
    for n = 1:Ntot
        spikediff = diff(spikes(spikes(:,2)==n,1))<burstISI;
        nsingle(n) = max(sum(spikediff == 0)+1,0);
        isidiff = diff(spikediff);
        nburst(n) = sum(isidiff==1);
        nsingle(n) = nsingle(n) - sum(isidiff==-1);
    end
    if isempty(data_Ne); pBursts_Ne(l1,l2) = NaN; else pBursts_Ne(l1,l2) = sum(nburst(1:Ne))/(sum(nsingle(1:Ne)) + sum(nburst(1:Ne))); end
    if isempty(data_Ni1); pBursts_Ni1(l1,l2) = NaN; else pBursts_Ni1(l1,l2) = sum(nburst(Ne+1:Ne+Ni1))/(sum(nsingle(Ne+1:Ne+Ni1)) + sum(nburst(Ne+1:Ne+Ni1))); end
    if isempty(data_Ni2); pBursts_Ni2(l1,l2) = NaN; else pBursts_Ni2(l1,l2) = sum(nburst(Ne+Ni1+1:Ntot))/(sum(nsingle(Ne+Ni1+1:Ntot)) + sum(nburst(Ne+Ni1+1:Ntot))); end
    
else
    
    tbins = tsel:dthist:Nt;
    LFP_sel = zeros(1,length(tt));
    
end