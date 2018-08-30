function pathTemp = remove_sys_files(pathTemp)

blCell = 1;
if ischar(pathTemp)
    pathTemp = {pathTemp};
    blCell = 0;
end

for zz = numel(pathTemp) : -1 : 1
    [~, fileTemp, ~] = fileparts(pathTemp{zz});
    if strcmpi(fileTemp(1:2), '._')
        pathTemp(zz) = [];
    end
end
clear zz

if blCell == 0
    if iscell(pathTemp)
       pathTemp = pathTemp{1}; 
    end
end