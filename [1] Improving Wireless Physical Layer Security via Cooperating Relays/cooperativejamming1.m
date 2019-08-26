%clear all; close all;
avgRsCJ=zeros(1,61);
c=3.5; %path loss exponent
%change Source-Eavesdropper distance
countNaN=0;
for dSE=30:90
    sumRsCJ=0;
    count=0;
  for itnr=1:100000
    e=(randn(3,1) + randn(3,1)*1i)/sqrt(2);
    PS=10^0.3;
    phi2=10^-6;
    %Source Destination
    dSD=50;
    hSD=((dSD)^(-c/2))*e(1,1);
    %Relays Destination
    dRD=25;
    hRD=((dRD)^(-c/2))*e;
    RRD=hRD*ctranspose(hRD);
   %Source Eavesdroppers
   hSE=((dSE)^(-c/2))*e(1,1);
   %Relays Eavesdroppers
   dRE=dSE-25;
   hRE=((dRE)^(-c/2))*e;
   RRE=hRE*ctranspose(hRE);
   %Optimal weight
   u=sqrt(PS)/norm((norm(hRD)^2)*hRE-ctranspose(hRD)*hRE*hRD);
   w=u*(norm(hRD)^2)*hRE-u*ctranspose(hRD)*hRE*hRD;
   %power of signal
   PRD=ctranspose(w)*RRD*w;
   PRE=ctranspose(w)*RRE*w;
   %secrecy rate
   RsCJ=log2(1+((PS*(norm(hSD)^2))/(phi2+PRD)))-log2(1+((PS*(norm(hSE)^2))/(phi2+PRE)));

   if (isnan(RsCJ)==1)
      countNaN=countNaN+1;
   elseif(RsCJ>=0)
      sumRsCJ=sumRsCJ+RsCJ;
      count=count+1;
   else
       sumRsCJ=sumRsCJ+RsCJ*0;
       count=count+1;    
   end
 end
 avgRsCJ(dSE-29)=sumRsCJ/count;
 hold on;
end
hold off;
dSE=30:90;
figure;
x=plot(dSE,avgRsCJ,'-x');
x.LineWidth=2;
x.Color = [0 0 0];
title('\fontsize{18}Average secrecy rate vs distance of source-eavesdropper channel');
xlabel('\fontsize{18}source-eavesdropper distance (m)');
ylabel('\fontsize{18}Secrecy rate (bits/s/Hz)');
disp(countNaN);

