function J=matrixJ(Ps,Rffavg,Pjavg,Rhhavg,phi2,I,Rggavg,PR);
J=inv(Ps*Rffavg+Pjavg*Rhhavg+phi2*I+PR*Rggavg);
end