clear
clc

fileID = fopen("..\..\gmsh-files\quarter-plate-with-hole-quad.msh", 'r');
if fileID == -1
    error('无法打开文件');
end


x_coor = [];
y_coor = [];
maxCols = 0; 
inDataSection = false;
 
tline = fgetl(fileID);

while ischar(tline)
    if strcmp(tline, '$Nodes')
        inDataSection = true;
        tline = fgetl(fileID);
        continue;
    elseif strcmp(tline, '$EndNodes')
        inDataSection = false;
    elseif inDataSection
        rowStrs = strsplit(tline);
        rowData = str2double(rowStrs);
        numRows = numel(rowData);
        if numRows > maxCols
            maxCols = numRows;
        end
        rowDataPadded = [rowData; zeros(maxCols - numRows, 1)]';
        rowDataPadded = rowDataPadded(1:maxCols);
        if numRows == 4
            x_coor = [x_coor; rowDataPadded(end-2)'];
            y_coor = [y_coor; rowDataPadded(end-1)'];
        end
    end
    tline = fgetl(fileID);
end


fileID = fopen("..\..\gmsh-files\quarter-plate-with-hole-quad.msh", 'r');

IEN = [];
maxCols = 0; 
inDataSection = false;

tline = fgetl(fileID);

while ischar(tline)
    if strcmp(tline, '$Elements')
        inDataSection = true;
        tline = fgetl(fileID);
        continue;
    elseif strcmp(tline, '$EndElements')
        inDataSection = false;
    elseif inDataSection
        rowStrs = strsplit(tline);
        rowData = str2double(rowStrs);
        numRows = numel(rowData);
        if numRows > maxCols
            maxCols = numRows;
        end
        rowDataPadded = [rowData; zeros(maxCols - numRows, 1)]';
        rowDataPadded = rowDataPadded(1:maxCols);
        if numRows == 9
            IEN = [IEN; rowDataPadded(end-3:end)'];
        end
    end
    tline = fgetl(fileID);
end


fclose(fileID);
