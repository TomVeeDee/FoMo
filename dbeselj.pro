
;Bessel function of the 1st kind (J)

FUNCTION dbeselj,X,n
  if n eq 0 then return, beselj(X,1,/double) else $
     return, -beselj(X,n+1,/double)+float(n)/X*beselj(X,n,/double)
end
