function hderi=hderivative(PR,Ps,hmatrix,JPs,Rffavg);
hderi=PR*(ctranspose(hmatrix)*JPs*hmatrix-Ps*ctranspose(hmatrix)*JPs^2*Rffavg*hmatrix);
end