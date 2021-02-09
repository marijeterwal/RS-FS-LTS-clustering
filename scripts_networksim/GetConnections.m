
%{

--- GetConnections ---
Defines connection matrix defined by the settings in GetSettings.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

% synapses
Syn     = zeros(Ne+Ni); %square matrix of synapses

nci     = [Ne Ni1 Ni2];
Nct     = length(nci);
cci     = cumsum([0 nci]);
indx    = {};
for i = 1:Nct
    indx{i} = 1:nci(i);
end

rng('default');
rng(seednr+4);

for i = 1:Nct %id receiver
    for j = 1:Nct %id sender
        Syn(cci(i)+indx{i},cci(j)+indx{j}) = SynConMat(i,j)*(rand(nci(i),nci(j))<pconSyn(i,j));
    end
end

Syn = sparse(Syn .* (AvgWeights + StdWeights*randn(Ntot, Ntot)));