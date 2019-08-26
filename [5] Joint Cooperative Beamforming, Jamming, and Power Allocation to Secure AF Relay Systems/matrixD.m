function D=matrixD(I,PR,APs,Rggavg)
D=I+PR*ctranspose(APs)*Rggavg*inv(APs);
end