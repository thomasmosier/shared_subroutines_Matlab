function pathOut = subset(pathIn, indRow, indCol, nm)

filesTemp = dir(fullfile(pathIn, '*.*'));
    filesTemp = struct2cell(filesTemp);
    filesTemp = filesTemp(1,end);
    
[~,~, ext] = fileparts(char(filesTemp{1}));
files = find_files(pathIn,ext(2:end));

mkdir(pathIn, nm);
pathOut = fullfile(pathIn, nm);

sMeta.currTime = [NaN, NaN];


for ii = 1 : numel(files(:))
	sData = read_geodata(files{ii}, sMeta, 'none', 'no_disp');
     
    [~,fileNm,ext]= fileparts(files{ii});
    
	sData.data = sData.data(indRow(1):indRow(end),indCol(1):indCol(end));
	sData.lat = sData.lat(indRow(1):indRow(end));
	sData.lon = sData.lon(indCol(1):indCol(end));
     
    indPrec = find(~isnan(sData.data) == 1 & sData.data ~= 0, 1, 'first');
    prec = precision(sData.data(indPrec));
    
    write_geodata(fullfile(pathOut,[nm '_' fileNm, ext]), sData, sMeta, prec, ext, 'no_disp');
end
