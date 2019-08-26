eavesdropper=25;
destination=50;
Ptdb=10;%Watt=db+30
Psdb=3;
Pt=10^(Ptdb/10);
Ps=10^(Psdb/10);
c=3.5/2;
phi2=10^-6; %Noise power
countNaN=0;
avgRs=zeros(1,51);
count0=0;
checkmax=0;
maxavgRs=0;
Pjmax=0;
Prmax=0;
%---------------------------------------------------------------------------
%random 8 intermediate nodes between 0 to 50m
a=1;
b=49;
%relay=round((b-a).*rand(8,1) + a) ;%random matrix [8 column, 1 row] for numbers between a and b
relay= sort(round((b-a)*rand(8,1)+a),'ascend')
e=rand + rand*1i;
%%%%%%%relay=transpose([2;11;11;12;19;43;43;46])
for Pjdb=-15:0.5:10
    Pj=10^(Pjdb/10);
    %=Ptdb-Pjdb-Psdb;
    Pr=Pt-Ps-Pj;
    %Pr=10^(Prdb/10)
    totalRs=0;
for time=1:100 %Check multiple times
%relay=round((b-a).*rand(8,1) + a) ;
hSE=(eavesdropper^(-c))*e;
%----------------------------------------------------------------------------------
    %Find best pair of relay m,n(m<n)
    for m=1:8
        for n=1:8
            if m~=n
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
                %thetaD2=(alpha^2)*(Pr/phi2)*(norm(hREn*hRDm-hREm*hRDn)^2);
                thetaD2=(Pr*norm(w*transpose(hRD)))/(phi2);
                
                if checkmax < thetaD2
                    checkmax=thetaD2;
                    mn=[m n];
                    hRED=[hREm hREn hRDm hRDn];
                end
            end
        end
    end
    %----------------------------------------------------------------------
    %final result of best relay m,n
    thetaD2=checkmax;
    checkmax=0;
    m=mn(1,1);
    n=mn(1,2);
    hREm=hRED(1,1);
    hREn=hRED(1,2);
    hRE=transpose([hREm hREn]);
    hRDm=hRED(1,3);
    hRDn=hRED(1,4);
    hRE=transpose([hRDm hRDn]);
%---------------------------------------------------------------------------
%Find best group of jammers o,p,q (o<p<q)

    for o=1:8
        if o~=m && o~=n
            for p=1:8
                if p~=m && p~=n && p~=o
                    for q=1:8
                        if q~=m && q~=n && q~=p && q~=o
                            %
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
                            hJE=transpose([hJEo hJEp hJEq]);
                                
                            if checkmax < (norm((u)*hJE))^2
                                checkmax=(norm((u)*hJE))^2;
                                opq=[o p q];
                            end
                        end
                    end
                end
            end
        end
    end
    %-----------------------------------------------------------------------
    %final result of jammers o,p,q
    thetaE1=(Ps*(norm(hSE)^2))/(phi2+Pj*checkmax);
    o=opq(1,1);
    p=opq(1,2);
    q=opq(1,3);
    Rs=1/2*log2(((1+thetaD2)/(1+thetaE1)));
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
avgRs(Pjdb*2+31)=totalRs/(time-countNaN);
if (maxavgRs<avgRs(Pjdb*2+31))
    maxavgRs=avgRs(Pjdb*2+31);
    Pjmax=Pjdb;
    Prmax=Ptdb-Psdb-Pjmax;
end
countNaN=0;
end
%---------------------------------------------------------------------------
%Plot result Rs vs Power of Jammers
Pjdb=-15:0.5:10;
%figure;
x=plot(Pjdb,avgRs,'-x');
x.LineWidth=2;
x.Color = [0 0 0];
title('\fontsize{16}Average secrecy rate vs power of jammers');
xlabel('\fontsize{16}Power of jammers (dB)');
ylabel('\fontsize{16}Secrecy rate (bits/s/Hz)');
disp(maxavgRs);
disp(Pjmax);
disp(Prmax);
disp(count0);