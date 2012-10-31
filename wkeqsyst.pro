FUNCTION wkeqsyst,X

; Dispersion relation in Edwin & Roberts 1983
; returns L(w,k) where L is dispersion relation = 0

  common vars1, A1, A2, A3, A4, B1, B2, B3, B4, n, ka
  dimx = n_elements(X)
  bslij = dblarr(dimx)
  bslky = dblarr(dimx)
  dbslij = dblarr(dimx)
  dbslky = dblarr(dimx)
  fracij = dblarr(dimx)
  fracky = dblarr(dimx)
  sigi = fltarr(dimx)
  sigk = fltarr(dimx)
  L = dblarr(dimx)

  nummoa2 = (ka^2*A1^2-X^2)*(ka^2*A2^2-X^2)
  dnummoa2 = (ka^2*A3^2-X^2)*(A1^2+A2^2)
  moa2 = nummoa2/dnummoa2
  nummea2 = (ka^2*B1^2-X^2)*(ka^2*B2^2-X^2)
  dnummea2 = (ka^2*B3^2-X^2)*(B1^2+B2^2)
  mea2 = nummea2/dnummea2

  reg1 = where(moa2 ge 0. and mea2 gt 0.)
  reg2 = where(moa2 ge 0. and mea2 lt 0.)
  reg3 = where(moa2 lt 0. and mea2 ge 0.)
  reg4 = where(moa2 lt 0. and mea2 le 0.)

  if reg1[0] ne -1 then begin
     bslij[reg1] = beseli(sqrt(abs(moa2[reg1])),n,/double)
     bslky[reg1] = beselk(sqrt(abs(mea2[reg1])),n,/double)
     dbslij[reg1] = dbeseli(sqrt(abs(moa2[reg1])),n)
     dbslky[reg1] = dbeselk(sqrt(abs(mea2[reg1])),n)
     fracij[reg1] = dbslij[reg1]/bslij[reg1]
     fracky[reg1] = dbslky[reg1]/bslky[reg1]
     sigi[reg1] = 1.
     sigk[reg1] = 1.
  endif
  if reg2[0] ne -1 then begin
     bslij[reg2] = beseli(sqrt(abs(moa2[reg2])),n,/double)
     bslky[reg2] = besely(sqrt(abs(mea2[reg2])),n,/double)
     dbslij[reg2] = dbeseli(sqrt(abs(moa2[reg2])),n)
     dbslky[reg2] = dbesely(sqrt(abs(mea2[reg2])),n)
     fracij[reg2] = dbslij[reg2]/bslij[reg2]
     fracky[reg2] = dbslky[reg2]/bslky[reg2]
     sigi[reg2] = 1.
     sigk[reg2] = -1.
  endif
  if reg3[0] ne -1 then begin
     bslij[reg3] = beselj(sqrt(abs(moa2[reg3])),n,/double)
     bslky[reg3] = beselk(sqrt(abs(mea2[reg3])),n,/double)
     dbslij[reg3] = dbeselj(sqrt(abs(moa2[reg3])),n)
     dbslky[reg3] = dbeselk(sqrt(abs(mea2[reg3])),n)
     fracij[reg3] = dbslij[reg3]/bslij[reg3]
     fracky[reg3] = dbslky[reg3]/bslky[reg3]
     sigi[reg3] = -1.
     sigk[reg3] = 1.
  endif
  if reg4[0] ne -1 then begin
     bslij[reg4] = beselj(sqrt(abs(moa2[reg4])),n,/double)
     bslky[reg4] = besely(sqrt(abs(mea2[reg4])),n,/double)
     dbslij[reg4] = dbeselj(sqrt(abs(moa2[reg4])),n)
     dbslky[reg4] = dbesely(sqrt(abs(mea2[reg4])),n)
     fracij[reg4] = dbslij[reg4]/bslij[reg4]
     fracky[reg4] = dbslky[reg4]/bslky[reg4]
     sigi[reg4] = -1.
     sigk[reg4] = -1.
  endif

  if reg1[0] ne -1 then L[reg1] = X[reg1]^2*(B4*sigi[reg1]*sqrt(abs(moa2[reg1]))*dbslij[reg1]*bslky[reg1]-A4*sigk[reg1]*sqrt(abs(mea2[reg1]))*dbslky[reg1]*bslij[reg1])-$
    ka^2*(B4*B2^2*sigi[reg1]*sqrt(abs(moa2[reg1]))*dbslij[reg1]*bslky[reg1]-A4*A2^2*sigk[reg1]*sqrt(abs(mea2[reg1]))*dbslky[reg1]*bslij[reg1])
  if reg2[0] ne -1 then L[reg2] = X[reg2]^2*(B4*sigi[reg2]*sqrt(abs(moa2[reg2]))*dbslij[reg2]*bslky[reg2]-A4*sigk[reg2]*sqrt(abs(mea2[reg2]))*dbslky[reg2]*bslij[reg2])-$
    ka^2*(B4*B2^2*sigi[reg2]*sqrt(abs(moa2[reg2]))*dbslij[reg2]*bslky[reg2]-A4*A2^2*sigk[reg2]*sqrt(abs(mea2[reg2]))*dbslky[reg2]*bslij[reg2])
  if reg3[0] ne -1 then L[reg3] = X[reg3]^2*(B4*sigi[reg3]*sqrt(abs(moa2[reg3]))*dbslij[reg3]*bslky[reg3]-A4*sigk[reg3]*sqrt(abs(mea2[reg3]))*dbslky[reg3]*bslij[reg3])-$
    ka^2*(B4*B2^2*sigi[reg3]*sqrt(abs(moa2[reg3]))*dbslij[reg3]*bslky[reg3]-A4*A2^2*sigk[reg3]*sqrt(abs(mea2[reg3]))*dbslky[reg3]*bslij[reg3])
  if reg4[0] ne -1 then L[reg4] = X[reg4]^2*(B4*sigi[reg4]*sqrt(abs(moa2[reg4]))*dbslij[reg4]*bslky[reg4]-A4*sigk[reg4]*sqrt(abs(mea2[reg4]))*dbslky[reg4]*bslij[reg4])-$
    ka^2*(B4*B2^2*sigi[reg4]*sqrt(abs(moa2[reg4]))*dbslij[reg4]*bslky[reg4]-A4*A2^2*sigk[reg4]*sqrt(abs(mea2[reg4]))*dbslky[reg4]*bslij[reg4])

  RETURN,L

;  RETURN,X^2*(B4*sqrt(moa2)*dbslij*bslky-A4*sqrt(mea2)*dbslky*bslij)-ka^2*(B4*B2^2*sqrt(moa2)*dbslij*bslky-A4*A2^2*sqrt(mea2)*dbslky*bslij)

;   RETURN,[X[3],$
;           X[0]^2*(A1^2+A2^2)*(X[3]^2*A3^2-X[2]^2)-(X[3]^2*A1^2-X[2]^2)*(X[3]^2*A2^2-X[2]^2),$
;           X[1]^2*(B1^2+B2^2)*(X[3]^2*B3^2-X[2]^2)-(X[3]^2*B1^2-X[2]^2)*(X[3]^2*B2^2-X[2]^2),$
;           X[2]^2*(B4*X[0]*dbslij*bslky-A4*X[1]*dbslky*bslij)-X[3]^2*(B4*B2^2*X[0]*dbslij*bslky-A4*A2^2*X[1]*dbslky*bslij)]
end
