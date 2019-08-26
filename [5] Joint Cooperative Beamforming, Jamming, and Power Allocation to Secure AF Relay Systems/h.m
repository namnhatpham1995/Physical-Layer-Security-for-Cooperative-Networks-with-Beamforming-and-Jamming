function hcoeff=h(PR,Ps,hmatrix,JPs);
hcoeff=PR*Ps*ctranspose(hmatrix)*JPs*hmatrix;
end