%{

--- ChangeParameters ---

Marije ter Wal - 2021
m.j.terwal@bham.ac.uk

%}

%  change loop variables to their current values
if ~isnan(duml1); eval([genvarname(loop1) '= duml1(l1);']); end
if ~isnan(duml2); eval([genvarname(loop2) '= duml2(l2);']); end
