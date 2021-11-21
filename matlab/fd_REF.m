fileID = fopen('fermi-dirac.txt','wt');
for lg2eta=-20:20
    fd=F(0,2^lg2eta,256);
    %disp(fd);
    %fprintf('%f',fd);
    fprintf(fileID,'%d\t%.16e\n',lg2eta,fd);
end
fclose(fileID);
winopen('fermi-dirac.txt');  % Only if using Windows!