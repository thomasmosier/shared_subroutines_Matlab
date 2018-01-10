function output = extract_field(S, name)
%Reproduces MAtlab function "extractfield" in mapping toolbox
nEntry = numel(S);
output = cell(nEntry, 1);
for ii = 1 : nEntry
    output{ii} = S(ii).(name);
end

