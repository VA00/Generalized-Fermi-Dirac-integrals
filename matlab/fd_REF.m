fileID = fopen('fermi-dirac.txt','wt');
k=0;
theta=4;
for lg2eta=-2:2:26 % from -2 step 2 to 26
    fd00=   F(k,2^lg2eta,theta);
    fd01=dF01(k,2^lg2eta,theta);
    fd02=dF02(k,2^lg2eta,theta);
    fd03=dF03(k,2^lg2eta,theta);
    fd10=dF10(k,2^lg2eta,theta);
    fd11=dF11(k,2^lg2eta,theta);
    fd12=dF12(k,2^lg2eta,theta);
    fd20=dF20(k,2^lg2eta,theta);
    fd21=dF21(k,2^lg2eta,theta);
    fd30=dF30(k,2^lg2eta,theta);
    fprintf(fileID,'%d\t%.17e\t%.17e\t%.17e\t%.17e\t',lg2eta,fd00,fd01,fd02,fd03);
    fprintf(fileID,'%.17e\t%.17e\t%.17e\t',fd10,fd11,fd12);
    fprintf(fileID,'%.17e\t%.17e\t',fd20,fd21);
    fprintf(fileID,'%.17e\n',fd30);
end
fclose(fileID);
winopen('fermi-dirac.txt');  % Only if using Windows!