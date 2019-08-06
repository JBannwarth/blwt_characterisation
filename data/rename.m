folderIn = '.';
files = dir( folderIn );

fileNames = {files.name}';
fileNames = fileNames( endsWith( fileNames, '.thA' ) );
newNames = strcat( erase( fileNames, ' (Ve).thA' ), '.thA' );
newNames = replace( newNames, 'course', 'coarse' );

for i = 1:length(newNames)
    movefile( fileNames{i}, newNames{i} )
end