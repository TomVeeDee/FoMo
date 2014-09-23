
; modified Bessel function of the second kind (K):

FUNCTION dbeselk,X,n
  if n eq 0 then return, -beselk(X,1,/double) else $
     return, -beselk(X,n-1,/double)-float(n)/X*beselk(X,n,/double)
end
