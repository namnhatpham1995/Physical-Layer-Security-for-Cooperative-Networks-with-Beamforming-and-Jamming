eavesdropper=25;
destination=50;
phi2=10^-6;
c=3.5/2;
a=1;
b=49;
maxRafi=zeros(1,13);
maxRdfi=zeros(1,13);
relay= sort(round((b-a)*rand(8,1)+a),'ascend');
e=rand + rand*1i;
for Psdb=1:13
    Ps=10^(Psdb/10);
    for ri=1:8
        hid=(abs(relay(ri,1)-destination)^(-c))*e;
        hsi=(abs(relay(ri,1))^(-c))*e;
        hie=(abs(relay(ri,1)-eavesdropper)^(-c))*e;
        Rafid=log2(1+(norm(hsi)^2*norm(hid)^2*Ps)/(2*(norm(hsi)^2+norm(hid)^2)*phi2));
        Rafie=log2(1+(norm(hsi)^2*norm(hie)^2*Ps)/(2*(norm(hsi)^2+norm(hie)^2)*phi2));
        
         Rafi=Rafid-Rafie;
         if (Rafi<0)
             Rafi=0;
         end
       
        
        if (maxRafi(1,Psdb)<Rafi)
            maxRafi(1,Psdb)=Rafi;
        end
        
    end
    for ri=1:8
        hid=(abs(relay(ri,1)-destination)^(-c))*e;
        hsi=(abs(relay(ri,1))^(-c))*e;
        hie=(abs(relay(ri,1)-eavesdropper)^(-c))*e;
        if norm(hsi)^2<norm(hid)^2;
            minhsihid=norm(hsi)^2;
        else
            minhsihid=norm(hid)^2;
        end
        Rdfsid=log2(1+(minhsihid*Ps)/(2*phi2));
        Rdfie=log2(1+(norm(hie)^2*Ps)/(2*phi2));
        Rdfi=Rdfsid-Rdfie;
         if (Rdfi<0)
             Rdfi=0;
         end
        if (maxRdfi(1,Psdb)<Rdfi)
            maxRdfi(1,Psdb)=Rdfi;
        end
    end
end
 Psdb=1:13;
 plot(Psdb,maxRafi,'-x');
 hold on;
  plot(Psdb,maxRdfi,'-o');
  hold off;
title('Average secrecy rate vs power of source');
xlabel('Power of source (dB)');
ylabel('Secrecy rate (bits/s/Hz)');
disp(maxRafi)
disp(maxRdfi)