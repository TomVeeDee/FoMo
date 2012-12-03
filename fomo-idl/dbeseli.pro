
;modified Bessel function of the 1st kind (I)

FUNCTION dbeseli,X,n
  if n eq 0 then return, beseli(X,1,/double) else $
     return, beseli(X,n+1,/double)+float(n)/X*beseli(X,n,/double)
end
