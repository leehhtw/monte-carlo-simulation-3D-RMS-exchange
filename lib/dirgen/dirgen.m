function gdir = dirgen(ndir)
filePath = matlab.desktop.editor.getActiveFilename;
root = fileparts(filePath);
system(sprintf('dirgen %u %s -cartesian -force',...
    ndir,...
    fullfile(root,'dirgen_temp.txt') ));
fid = fopen(fullfile(root,'dirgen_temp.txt'),'r');
fgetl(fid);
gdir = zeros(ndir,3);
for i = 1:ndir
    tline = fgetl(fid);
    gdir(i,:) = str2num(tline);
end
fclose(fid);
end