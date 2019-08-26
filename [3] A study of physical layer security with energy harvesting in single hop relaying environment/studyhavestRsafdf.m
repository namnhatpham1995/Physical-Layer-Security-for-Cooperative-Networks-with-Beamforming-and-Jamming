dSR=10;
dRD=15;
dBS=7;
dBR=7;
phi2=10^-6;
Pb=10^3;
alpha=0.99;
nguy=0.9;
e=(randn(1,1) + randn(1,1)*1i);
c=3.5;
Pj=10;
Rsaf=zeros(1,26);
Rsdf=zeros(1,26);
avgRsdf=zeros(1,26);
avgRsaf=zeros(1,26);
for time=1:10000
    e=(randn(1,1) + randn(1,1)*1i);
    for dRE=15:40

        dSE=sqrt(dSR^2+dRE^2);

        hSR=(dSR^(-c/2))*e;
        hRD=(dRD^(-c/2))*e;
        hBS=(dBS^(-c/2))*e;
        hBR=(dBR^(-c/2))*e;
        hRE=(dRE^(-c/2))*e;
        hSE=(dSE^(-c/2))*e;
        Ps=(2*nguy*Pb*(norm(conj(hBS))^2)*alpha)/(1-alpha);
        Pr=(2*nguy*Pb*(norm(conj(hBR))^2)*alpha)/(1-alpha);
        G=1/sqrt(Ps*(norm(hSR)^2));
        alphaRD=(norm(hRD)^2)/phi2;
        alphaSE=(norm(hSE)^2)/phi2;
        alphaRE=(norm(hRE)^2)/phi2;


        Rddf=(1/2)*log2(1+Pr*alphaRD);
        Redf=(1/2)*log2(1+((Ps*alphaSE)/(1+Pj*alphaRE))+((Pr*alphaRE)/(1+Pj*alphaSE)));
        Rsdf(1,dRE-14)=Rddf-Redf;
        if Rsdf(1,dRE-14)<0
            Rsdf(1,dRE-14)=0;
        end

        Rdaf=(1/2)*log2(1+G^2*Ps*alphaRD);
        Reaf=(1/2)*log2(1+((Ps*alphaSE)/(1+Pj*alphaRE))+((G^2*Ps*alphaRE)/(1+Pj*alphaSE)));
        Rsaf(1,dRE-14)=Rdaf-Reaf;
        if Rsaf(1,dRE-14)<0
            Rsaf(1,dRE-14)=0;
        end
        
        %avgRsdf(1,dRE-14)=avgRsdf(1,dRE-14)+Rsdf(1,dRE-14);
        %avgRsaf(1,dRE-14)=avgRsaf(1,dRE-14)+Rsaf(1,dRE-14);
    end
    avgRsdf=avgRsdf+Rsdf;
    avgRsaf=avgRsaf+Rsaf;
end
avgRsdf=avgRsdf/time;
avgRsaf=avgRsaf/time;
%for dRE=15:40
%avgRsaf(1,dRE-14)=avgRsaf(1,dRE-14)/time;
%avgRsdf(1,dRE-14)=avgRsdf(1,dRE-14)/time;
%end
dRE=15:40;
x=plot(dRE,avgRsdf,'-o');
x.LineWidth=2; 
hold on;
y=plot(dRE,avgRsaf,'-x');
y.LineWidth=2; 
hold off;
title('\fontsize{18}Average secrecy rate vs relay-eavesdropper distance');
xlabel('\fontsize{18}Relay-eavesdropper distance (m)');
ylabel('\fontsize{18}Secrecy rate (bits/s/Hz)');
