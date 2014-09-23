
;Bessel function of the 1st kind (J)

FUNCTION dbeselj,X,n  
  if n eq 0 then begin 
     return, -beselj(X,1,/double) 
  endif else begin
     if min(abs(X)) eq 0. then return, 0.5*(-beselj(X,n+1,/double)+beselj(X,n-1,/double)) else return, beselj(X,n-1,/double)-float(n)/X*beselj(X,n,/double)
  endelse
end
