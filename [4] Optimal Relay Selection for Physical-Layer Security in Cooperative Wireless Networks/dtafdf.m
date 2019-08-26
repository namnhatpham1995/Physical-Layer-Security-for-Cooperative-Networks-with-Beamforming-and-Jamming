alphasi=1;
alphaid=1;
alphaie=1;
%alphaid=sigid^2/sigsd^2;
%alphaie=sigie^2/sigse^2;
%MER
%lamdade=sigsd^2/sigse^2;

for MER=-5:15
    lamdade=10^(MER/10);
    Pdt(MER+6)=1/(1+lamdade);
    Paf(MER+6)=alphaie/(alphaie+alphaid*lamdade);
    Pdf(MER+6)=(alphaid+alphasi)/(alphaid+alphasi+alphasi*alphaid*alphaie^-1*lamdade);
    
end
figure;
MER=-5:15
x=semilogy(MER,Pdt,'-^');
x.LineWidth=2;
hold on;
xx=semilogy(MER,Paf.^2,'-x');
xx.LineWidth=2;
xxx=semilogy(MER,Pdf.^2,'-d');
xxx.LineWidth=2;
xxxx=semilogy(MER,Paf.^4,'-s');
xxxx.LineWidth=2;
xxxxx=semilogy(MER,Pdf.^4,'--');
xxxxx.LineWidth=2;
xxxxx.Color = [0 0 0];
hold off;
title('\fontsize{18}Direct transmission vs P-AFbORS vs P-DFbORS');
ylabel('\fontsize{18}Intercept probability');
xlabel('\fontsize{18}Main-to-eavesdropper ratio MER(dB)');