%prompt = 'How many relays N do you want? ';
%n = input(prompt) %numbers of relay N
n=8;
Psdb=0;
Ps=10^(Psdb/10);
phi2=10^-6;
fRn=zeros(1,n);
gRn=zeros(1,n);
hRn=zeros(1,n);
cEn=zeros(1,n);
acf=zeros(1,n);
acg=zeros(1,n);
afg=zeros(1,n);
agh=zeros(1,n);
Rgg=zeros(1,n);
Rff=zeros(1,n);
Rcc=zeros(1,n);
Rhh=zeros(1,n);
Rs=zeros(1,9);
%I=eye(n-2,n-2)
%I=eye(2,2)
I=eye(n,n)
for count=1:n
%fRn(1,count)=normrnd(0,1)+normrnd(0,1)*i ;%row vector
%gRn(1,count)=normrnd(0,1)+normrnd(0,1)*i ;%row vector
%hRn(1,count)=normrnd(0,1)+normrnd(0,1)*i ;%row vector
%cEn(1,count)=normrnd(0,1)+normrnd(0,1)*i ;%row vector
%fE=normrnd(0,1)+normrnd(0,1)*i;
%qE=normrnd(0,1)+normrnd(0,1)*i;
fRn(1,count)=rand+rand*i ;%row vector
gRn(1,count)=rand+rand*i ;%row vector
hRn(1,count)=rand+rand*i ;%row vector
cEn(1,count)=rand+rand*i ;%row vector
fE=randn+randn*i;
qE=randn+randn*i;
acf(1,count)=cEn(1,count)*fRn(1,count);
acg(1,count)=cEn(1,count)*gRn(1,count);
afg(1,count)=fRn(1,count)*gRn(1,count);
agh(1,count)=gRn(1,count)*hRn(1,count);
Rff(1,count)=norm(fRn(1,count))^2;
Rgg(1,count)=norm(gRn(1,count))^2;
Rcc(1,count)=norm(cEn(1,count))^2;
Rhh(1,count)=norm(hRn(1,count))^2;
end
Rff=diag(Rff)
Rgg=diag(Rgg)
Rhh=diag(Rhh)
Rcc=diag(Rcc)
fRn=transpose(fRn)
hRn=transpose(hRn)
cEn=transpose(cEn)
gRn=transpose(gRn)



acf=transpose(acf);
acg=transpose(acg);
afg=transpose(afg);
agh=transpose(agh);

Rfg=afg*ctranspose(afg)
Rcg=acg*ctranspose(acg);
Rcf=acf*ctranspose(acf);

%H=transpose([acf agh])
%Hkernel=transpose(null(H,'r'))%find kernel/null space of matrix
H=[acf agh];
%Hnull=(Hkernel*inv(transpose(Hkernel)*Hkernel)*transpose(Hkernel))*H%projection of H onto nullspace
%Hnull=null(ctranspose(H))
Hnull=(H*inv(transpose(H)*H)*transpose(H))
%a=ctranspose(H)*Hnull
Rfgavg=ctranspose(Hnull)*Rfg*(Hnull);
Rggavg=ctranspose(Hnull)*Rgg*(Hnull);
Rffavg=ctranspose(Hnull)*Rff*(Hnull)
Rhhavg=ctranspose(Hnull)*Rhh*(Hnull);

for Pm=8:4:40

    PRdb=Pm*(n-1)/n;
    PR=10^(PRdb/10);
    Pjavgdb=PRdb/(n-1);
    Pjavg=10^(Pjavgdb/10);
     hmatrix=ctranspose(Hnull)*afg;
    a=phi2*Pjavg*norm(qE)^2;
    b=norm(fE)^2;
    
    JPs=matrixJ(Ps,Rffavg,Pjavg,Rhhavg,phi2,I,Rggavg,PR);
    J0=matrixJ(0,Rffavg,Pjavg,Rhhavg,phi2,I,Rggavg,PR);    
    hPs=h(PR,Ps,hmatrix,JPs);
    h0=h(PR,0,hmatrix,J0);
   
    hPsderi=hderivative(PR,Ps,hmatrix,JPs,Rffavg);
    h0deri=hderivative(PR,0,hmatrix,J0,Rffavg);
    m0=matrixm(a,b,h0,phi2,0,h0deri);
    mPs=matrixm(a,b,hPs,phi2,Ps,hPsderi);
%    if (m0<=0)
 %       Psdb=0;
 %       Ps=10^(Psdb/10);
%    else%if (mPs>=0)
%        Psdb=30;
%        Ps=10^(Psdb/10);
%    else
%        solve(matrixJ(Ps,Rffavg,Pjavg,Rhhavg,phi2,I,Rggavg,PR)
%        solve(h(PR,Pc,hmatrix,JPc),Pc)
%        solve(matrixm(a,b,h0,phi2,Pc,h0deri),Pc);
 %   end
    
    TavgPs=matrixTavg(Ps,Rffavg,Pjavg,Rhhavg,phi2,I);
    rowTavg=size(TavgPs);
    for count=1:rowTavg(1,1)
        TavgPs(count,count)=real(TavgPs(count,count));
    end

    APs=matrixA(TavgPs);
    DPs=matrixD(I,PR,APs,Rggavg);
    alpha1=inv(DPs)*ctranspose(inv(APs));
    alpha=norm(inv(DPs)*ctranspose(inv(APs))*ctranspose(Hnull)*afg);

    vopti=alpha*sqrt(PR)*inv(APs)*inv(DPs)*inv(ctranspose(APs))*ctranspose(Hnull)*afg;
     wopti=Hnull*vopti;
     Rd=(1/2)*log2(1+(Ps*ctranspose(wopti)*Rfg*wopti)/(phi2*(1+ctranspose(wopti)*Rgg*wopti)));
    %Rd=(1/2)*log2(1+(Ps*ctranspose(afg)*wopti*ctranspose(wopti)*afg)/(phi2*(1+ctranspose(wopti)*Rgg*wopti)));
     Re=(1/2)*log2(1+(Ps*norm(fE)^2)/(phi2*Pjavg*norm(qE)^2));
     if (Rd-Re<0)
        Rs(1,Pm/4-1)=0;
     else
        Rs(1,Pm/4-1)=Rd-Re
     end
end
  Pm=8:4:40;
x=plot(Pm,(Rs),'-x');
x.LineWidth=2;
x.Color = [0 0 0];
title('\fontsize{18}Average secrecy rate vs power relay and jammer');
xlabel('\fontsize{18}Power relay and jammer Pm (dB)');
ylabel('\fontsize{18}Secrecy rate (bits/s/Hz)');

 %plot(Pm,Rs+1+abs(Rs(1,1)),'-x');