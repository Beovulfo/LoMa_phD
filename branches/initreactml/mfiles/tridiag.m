function A = tridiag(a,b,c,n)

diagonal = a*ones(1,n);
infdiag  = b*ones(1,n-1);
supdiag  = c*ones(1,n-1);

A = diag(diagonal,0) + diag(infdiag,-1) + diag(supdiag,1);

end