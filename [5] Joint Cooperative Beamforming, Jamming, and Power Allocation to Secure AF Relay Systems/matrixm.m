function m=matrixm(a,b,hPs,phi2,Ps,hPsderi);
m=(a+b*Ps)*hPsderi-b*phi2-b*hPs;
end