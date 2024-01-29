
%{

--- GetConnections ---
Defines connection matrix defined by the settings in GetSettings.

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

nregions = length(Ne);
SynStore = {};

% connections within each local circuit
for nr = 1:nregions

    % synapses
    Syndum     = zeros(Ne(nr)+Ni(nr)); %square matrix of synapses

    nci     = [Ne(nr) Ni1(nr) Ni2(nr)];
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
            Syndum(cci(i)+indx{i},cci(j)+indx{j}) = SynConMat{nr,nr}(i,j)*(rand(nci(i),nci(j))<pconSyn{nr,nr}(i,j));
        end
    end

    SynStore{nr,nr} = sparse(Syndum .* (AvgWeights + StdWeights*randn(Ntot(nr), Ntot(nr))));
end

% connections between circuits
for nr1 = 1:nregions
    for nr2 = 1:nregions

        if nr1==nr2
            continue
        end

        % synapses
        Syndum     = zeros(Ne(nr1)+Ni(nr1),Ne(nr2)+Ni(nr2)); %square matrix of synapses

        nci1     = [Ne(nr1) Ni1(nr1) Ni2(nr1)];
        Nct1     = length(nci1);
        cci1     = cumsum([0 nci1]);
        indx1    = {};
        for i = 1:Nct1
            indx1{i} = 1:nci1(i);
        end

        nci2     = [Ne(nr2) Ni1(nr2) Ni2(nr2)];
        Nct2     = length(nci2);
        cci2     = cumsum([0 nci2]);
        indx2    = {};
        for i = 1:Nct2
            indx2{i} = 1:nci2(i);
        end

        rng('default');
        rng(seednr+4);

        for i = 1:Nct1 %id receiver
            for j = 1:Nct2 %id sender
                Syndum(cci1(i)+indx1{i},cci2(j)+indx2{j}) = SynConMat{nr1,nr2}(i,j)*(rand(nci1(i),nci2(j))<pconSyn{nr1,nr2}(i,j));
            end
        end

        SynStore{nr1,nr2} = sparse(Syndum .* (AvgWeights + StdWeights*randn(Ntot(nr1), Ntot(nr2))));
    end
end

% merge connection matrix
Syn = cell2mat(SynStore);