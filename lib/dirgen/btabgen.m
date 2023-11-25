function btab = btabgen(gdir,bval,t,root)
Nt = size(t,1);
Nb = numel(bval);
Ng = size(gdir,1);
btab = [];
for i = 1:Nt
    ti = t(i,:);
    btab = cat(1,btab,[0 0 0 0 ti]);
    for j = 1:Nb
        bj = bval(j);
%         gvalj = sqrt(bj/ti(2)^2/(ti(1)-ti(2)/3));
        for k = 1:Ng
            gdirk = gdir(k,:);
            btab = cat(1,btab,[gdirk bj ti]);
        end
    end
end
if nargin>3
fid = fopen(fullfile(root),'w');
for i = 1:size(btab,1)
    fprintf(fid,'%.16f %.16f %.16f %.16f %f %f %f \n',btab(i,:));
end
fclose(fid);
end
end