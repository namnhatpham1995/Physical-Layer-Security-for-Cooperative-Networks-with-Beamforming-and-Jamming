%prompt = 'How many relays N do you want? ';
%n = input(prompt) %numbers of relay N
n=5
destination=50;
c=3.5/2;
a=1;
b=49;
%relay= sort(round((b-a)*rand(n,1)+a),'ascend');%8 11 21 29 37
relay=[3;30;33;35;37];%37(jammer) 33(os) 35(oscj) 30(osj)
%e=rand + rand*1i;
Rsosj=zeros(1,49);
Rsos=zeros(1,49);
Rsoscj=zeros(1,49);
totalosj=zeros(1,49);
totalos=zeros(1,49);
totaloscj=zeros(1,49);
Pdb=50;
P=10^(Pdb/10);
for times=1:10000
    e=rand + rand*1i;
    for eavesdropper=1:49
        for m=1:n
            for o=1:n
                if o~=m
                   hRD=(abs(relay(m,1)-destination)^(-c));
                   hRE=(abs(relay(m,1)-eavesdropper)^(-c));
                   hJD=(abs(relay(o,1)-destination)^(-c));
                   hJE=(abs(relay(o,1)-eavesdropper)^(-c));
                   thetaRD=P*(norm(hRD)^2);
                   thetaRE=P*(norm(hRE)^2);
                   thetaJD=P*(norm(hJD)^2);
                   thetaJE=P*(norm(hJE)^2);

                   Rdosj=(1/2)*log2(1+(thetaRD/(1+thetaJD)));  %OSJ
                   Reosj=(1/2)*log2(1+(thetaRE/(1+thetaJE)));
                   Rsosj1=Rdosj-Reosj;
                   if Rsosj(1,eavesdropper)<Rsosj1
                   Rsosj(1,eavesdropper)=Rsosj1;
                   chosenjammerosj=o;
                   chosenrelayosj=m;
                   end

                   if thetaJE>thetaJD                           %OS
                       %Rdos=(1/2)*log2(1+(thetaRD/(1+thetaJD)));  
                       %Reos=(1/2)*log2(1+(thetaRE/(1+thetaJE)));
                       Rdos=(1/2)*log2(1+(thetaRD/thetaJD));  
                       Reos=(1/2)*log2(1+(thetaRE/thetaJE));
                       Rsos1=Rdos-Reos;
                       if Rsos(1,eavesdropper)<Rsos1
                       Rsos(1,eavesdropper)=Rsos1;
                       chosenjammeros=o;
                       chosenrelayos=m;
                       end 
                   end

                   Rdoscj=(1/2)*log2(1+thetaRD);
                   Reoscj=(1/2)*log2(1+(thetaRE/(1+thetaJE)));
                   Rsoscj1=Rdoscj-Reoscj;
                   if Rsoscj(1,eavesdropper)<Rsoscj1
                   Rsoscj(1,eavesdropper)=Rsoscj1;
                   chosenjammeroscj=o;
                   chosenrelayoscj=m;
                   end
                end
            end
        end    
    end
    totaloscj=totaloscj+Rsoscj;
    totalos=totalos+Rsos;
    totalosj=totalosj+Rsosj;
end
totaloscj=totaloscj/times;
totalos=totalos/times;
totalosj=totalosj/times;
eavesdropper=1:49;
%figure;
x=plot(eavesdropper,totalosj,'-x');
x.LineWidth=2; 
hold on;
y=plot(eavesdropper,totalos,'-o');
y.LineWidth=2; 
z=plot(eavesdropper,totaloscj,'--');
z.LineWidth=2; 
z.Color = [0 0 0];  
title('\fontsize{18}Average secrecy rate vs eavesdropper postion');
xlabel('\fontsize{18}eavesdropper position');
ylabel('\fontsize{18}Secrecy rate (bits/s/Hz)');
hold off;
disp('end')