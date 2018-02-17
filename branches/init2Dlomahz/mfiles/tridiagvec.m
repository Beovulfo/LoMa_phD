function A = tridiagvec(diagn,infdiag,supdiag)

N = length(diagn);

A = diag(diagn) + diag(infdiag,-1) + diag(supdiag,1);

end