destination=50;
Ptdb=40;%Watt=db+30
Psdb=34;
Pt=10^(Ptdb/10);
Ps=10^(Psdb/10);
c=3.5;
phi2=10^-6; %Noise power
countNaN=0;
avgRs=zeros(1,45);
count0=0;
%m= ;
%n= ;
%o= ;
%p= ;
%q= ;
%---------------------------------------------------------------------------
%random 8 intermediate nodes between 0 to 50m
%a=0;
%b=50;
%relay=round((b-a).*rand(8,1) + a) ;%random matrix [8 column, 1 row] for numbers between a and b
%relay= sort(round((b-a)*rand(8,1)+a),'ascend')
%e=rand + rand*1i
%e=rand
%relay=transpose([15 16 17 23 30 31 32 44])
%relay=transpose([1 4 13 14 15 21 25 40])
%%%relay=transpose([1 5 46 15 22 31 32 49])
%%%relay=transpose([2 5 13 18 40 44 45 49])
%relay=transpose([18 19 22 28 30 30 34 34])


for eavesdropper=5:49
    Pjdb=3;
    Pj=10^(Pjdb/10);
    Prdb=Ptdb-Psdb-Pjdb;
    Pr=10^(Prdb/10);
    totalRs=0;
for time=1:1000 %Check multiple times
hSE=(eavesdropper^(-c))*e;
%----------------------------------------------------------------------------------
hREm=(abs(relay(m,1)-eavesdropper)^(-c))*e;%dREm
hREn=(abs(relay(n,1)-eavesdropper)^(-c))*e;%dREn
hRDm=(abs(relay(m,1)-destination)^(-c))*e;%dRDm
hRDn=(abs(relay(n,1)-destination)^(-c))*e;%dRDn
hRD=[hRDm hRDn];
                %
alpha=1/((norm(hREm)^2)+(norm(hREn)^2));
wm=alpha*hREn;
wn=-alpha*hREm;
w=[wm wn];
                %
thetaD2=(Pr*norm(w*transpose(hRD)))/(phi2);
%---------------------------------------------------------------------------
hRpRm=(abs(relay(m,1)-relay(p,1))^(-c))*e;%(dRpRm)
hRpRn=(abs(relay(n,1)-relay(p,1))^(-c))*e;%(dRpRn)
hRqRn=(abs(relay(n,1)-relay(q,1))^(-c))*e;%(dRqRn)
hRqRm=(abs(relay(m,1)-relay(q,1))^(-c))*e;%(dRqRm)
hRoRm=(abs(relay(m,1)-relay(o,1))^(-c))*e;%(dRoRm)
hRoRn=(abs(relay(n,1)-relay(o,1))^(-c))*e;%(dRoRn)
u1=(hRpRm*hRqRn-hRpRn*hRqRm)/sqrt((hRpRm*hRqRn-hRpRn*hRqRm)^2+(hRoRn*hRqRm-hRoRm*hRqRn)^2+(hRoRm*hRpRn-hRpRm*hRoRn)^2);
u2=(hRoRn*hRqRm-hRoRm*hRqRn)/sqrt((hRpRm*hRqRn-hRpRn*hRqRm)^2+(hRoRn*hRqRm-hRoRm*hRqRn)^2+(hRoRm*hRpRn-hRpRm*hRoRn)^2);
u3=(hRoRm*hRpRn-hRpRm*hRoRn)/sqrt((hRpRm*hRqRn-hRpRn*hRqRm)^2+(hRoRn*hRqRm-hRoRm*hRqRn)^2+(hRoRm*hRpRn-hRpRm*hRoRn)^2);
u=([u1 u2 u3]);
hJEo=(abs(relay(o,1)-eavesdropper)^(-c))*e;%dJEo
hJEp=(abs(relay(p,1)-eavesdropper)^(-c))*e;%dJEp
hJEq=(abs(relay(q,1)-eavesdropper)^(-c))*e;%dJEq
hJE=([hJEo hJEp hJEq]);

thetaE1=(Ps*(norm(hSE)^2))/(phi2+Pj*checkmax);

Rs=(1/2)*log2((1+thetaD2)/(1+thetaE1));
%-----------------------------------------------------------------------------------------------
    %summarize Rs
    if(Rs>=0)
      totalRs = totalRs+Rs;
    else
        countNaN=countNaN+1;
        count0=count0+1;
    end
    
end
%---------------------------------------------------------------------------
%Average Rs at appropriate Pjdb
avgRs(eavesdropper-4)=totalRs/(time-countNaN);
countNaN=0;
end
%---------------------------------------------------------------------------
%Plot result Rs vs Power of Jammers
eavesdropper=5:49;
figure;
x=plot(eavesdropper,avgRs,'-x');
x.LineWidth=2;
x.Color = [0 0 0];
title('\fontsize{16}Average secrecy rate vs eavesdropper position');
xlabel('\fontsize{16}eavesdropper position');
ylabel('\fontsize{16}Secrecy rate (bits/s/Hz)');
disp(count0);