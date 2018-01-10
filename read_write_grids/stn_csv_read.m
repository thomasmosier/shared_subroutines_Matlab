function [dateIn, dataIn] = stn_csv_read(pathIn, stn)

%Read header:
fid = fopen(pathIn,'r');
hdrTemp = textscan(fid, '%[^\n]', 1);
fclose(fid);
hdrTemp = char(hdrTemp{1});

%Parse header:
indCom = regexpi(hdrTemp, ',');
if indCom(end) == numel(hdrTemp(:))
    nCol = numel(indCom);
else
    nCol = numel(indCom)+1;
end

hdr = cell(1, nCol);
for ii = 1 : numel(indCom)
    if ii == 1
        hdr{ii} = strrep(hdrTemp(1:indCom(ii)-1), ' ', '');
    else
        hdr{ii} = strrep(hdrTemp(indCom(ii-1)+1:indCom(ii)-1), ' ', '');
    end
end
if numel(hdr) > numel(indCom)
    hdr{end} = strrep(hdrTemp(indCom(end)+1:end), ' ', '');
end

indStn = find(strcmpi(hdr, stn) == 1);
    
if isempty(indStn)
   error('cropWaterReq:noStnFoundHist', ['Station ' stn ' not found in historical data file.']); 
end
    
dataArrayTemp = csvread(pathIn, 1, 0);


if strcmpi(hdr{1} , 'year') && strcmpi(hdr{2} , 'month') && strcmpi(hdr{3} , 'day') 
    dateIn = dataArrayTemp(:,1:3);
elseif strcmpi(hdr{1} , 'year') && strcmpi(hdr{2} , 'month') 
    dateIn = dataArrayTemp(:,1:2);
else
    error('cropWaterReq:noStnFoundHist', 'date representation not recognized');
end

dataIn = dataArrayTemp(:,indStn);