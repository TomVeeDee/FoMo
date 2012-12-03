
;Bessel function of the 2nd kind (Y)

FUNCTION dbesely,X,n
  if n eq 0 then return, besely(X,1,/double) else $
     return, -besely(X,n+1,/double)+float(n)/X*besely(X,n,/double)
end
