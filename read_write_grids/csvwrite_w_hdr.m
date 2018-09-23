function csvwrite_w_hdr(pathIn,M,hdr,varargin)

row = 0;
col = 0;
if ~isempty(varargin(:))
    if numel(varargin(:)) == 2
        row = varargin{1};
        col = varargin{2};
    else
       error('csvwriteWHdr:wrongVariableInput', ['The variable input argument has ' num2str(numel(varargin(:))) ' entries. Either 0 or 2 are expected.']) 
    end
end

if ~iscell(hdr)
    error('csvwriteWHdr:wrongHdrFormat', ['The header is class ' class(hdr) ', but cell array is expected.']);
end
textHeader = strjoin(hdr, ','); %insert commaas
textHeader = textHeader(:)';

[fold, file, ext] = fileparts(pathIn);
if isempty(ext)
    pathIn = fullfile(fold, [file, '.csv']);
end

%write header to file
fid = fopen(pathIn,'w+'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);

%Write data to file
if ismatrix(M) && ~iscell(M)
    dlmwrite(pathIn, M(row+1:end, col+1:end),'-append');
elseif iscell(M)
    fid = fopen(pathIn,'w'); 

    for ii = 1 : numel(M(:,1))
        dataTmp = cell(numel(M(1,:)), 1);
        for jj = 1 : numel(M(1,:))
            if isnumeric(M{ii,jj})
                dataTmp{jj} = num2str(M{ii,jj});
            else
                dataTmp{jj} = M{ii,jj};
            end
        end
        
        fprintf(fid, '%s\n', strjoin(dataTmp, ','));
    end
    fclose(fid);
else
    error('csvwriteWHdr:wrongDataFormat', ['The data is class ' class(M) ', but matrix or cell array is expected.']);
end
    