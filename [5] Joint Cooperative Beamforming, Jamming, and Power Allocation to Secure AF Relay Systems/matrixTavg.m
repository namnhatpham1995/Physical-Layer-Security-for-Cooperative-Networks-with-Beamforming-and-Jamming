function Tavg=matrixTavg(x,Rffavg,Pjavg,Rhhavg,phi2,I)
Tavg=x*Rffavg+Pjavg*Rhhavg+phi2*I;
end