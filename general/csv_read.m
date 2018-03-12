function [data, hdr] = csv_read(path)
%Read header and data from csv file

%Get header:
fileID = fopen(path);
tline = fgetl(fileID);
hdr = strtrim(strsplit(tline,','));
fclose(fileID);

%Get numbers:
data = csvread(path,1, 0);