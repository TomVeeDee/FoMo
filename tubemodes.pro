
PRO tubemodes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, waka_ini=waka_ini, waka_root, ka_root, normcase = normcase

if n_params(0) lt 1 then begin
   print,'tubemode, ro=ro, re=re, va=va, vae=vae, co=co, ce=ce, waka_ini=waka_ini, waka_root, ka_root, [normcase = normcase]'
   return
endif

; Solves dispersion relation based on parameters ro, re, va, vae, co,
; ce, waka_ini.
; INPUT:
; ro = internal density: rho[dimx/2,dimy/2,dimz/2]
; re = external density: rho[0,0,0]
; va = internal Alfven velocity: calculated at (dimx/2,dimy/2,dimz/2)
; vae = external Alfven velocity: calculated at (0,0,0)
; co = internal sound speed: calculated at (dimx/2,dimy/2,dimz/2)
; ce = external sound speed: calculated at (0,0,0)
; waka_ini = initial guess of w/k for solving dispersion relation
; OUTPUT:
; waka_root = roots of w/k from dispersion relation
; ka_root = values of ka corresponding to w/k roots (a = radius of cylinder)
; OPTIONAL:
; if keyword 'normcase' is set then dispersion relation is solved with
; adhoc values given below
; CALLS:
; function 'wkeqsyst.pro': dispersion relation

; define parameters:
common vars1, A1, A2, A3, A4, B1, B2, B3, B4, n, ka

n = 0 ; sausage mode

;xm = round(dimx/2.)
;ym = round(dimy/2.)
;zm = round(dimz/2.)

gamma = 5./3.
mup=1.25663706*1.e-6
nka = 201 ; number of roots
ka_ar = findgen(nka)/(nka-1.)*3.+1. ; initial values for ka
ka_root = fltarr(nka)

if keyword_set(normcase) then begin
   co = 1.
   ce = 0.5 * co
   va = 2 * co
   vae = 5 * co
   ro = 1.
   re = 0.1
endif

ct = co*va/sqrt(co^2+va^2)
cte = ce*vae/sqrt(ce^2+vae^2)

A1 = co & A2 = va & A3 = ct & A4 = ro 
B1 = ce & B2 = vae & B3 = cte & B4 = re 
mrk = 0
;step = waka_ini*0.005
waka_j = waka_ini

for j=0,nka-1 do begin
   ka_arrev = reverse(ka_ar)
;   ka = ka_ar[-1-j]
   ka = ka_arrev[j]
   wa0 = waka_j*ka
;   wa0 = (waka_ini+step*j)*ka
;   if mrk eq 1 then begin
;      root_1 = (wa_rooti((uniq(wa_rooti))[wrlocs])/ka)[0]
;      if root_1 gt (vae-0.02*vae) then step = waka_ini*0.0015
;      wa0 = (waka_ini+step*j)*ka
;   endif
   wa_rooti = 0.
   dwa1 = 1.
   k = 0
   while(dwa1 gt 1.e-4) do begin
      stepg = 0.001
      wa = wa0-stepg*wa0*k
      wa_2 = wa0-0.1*ka-stepg*wa0*k
      wa_3 = wa0-0.2*ka-stepg*wa0*k
; initial guess: iX[0]=moa, iX[1]=mea, iX[2]=wa, iX[3]=ka
      iX=[wa,wa_2,wa_3]
      wa1 = float(real_part(FX_ROOT(iX,'wkeqsyst',/double)))
;      wa1 = float(bisection(iX[0],'wkeqsyst',tol=1.e-4,radius=1.))
      if k eq 0 then wa_rooti = wa1 else wa_rooti = [wa_rooti,wa1]
      dwa1 = abs(va*ka-wa1)
      k = k+1
      if wa1/ka gt 0.98*vae then begin
         waka_j=va*1.1
         mrk = 1
      endif

;      if (wa1/ka lt 8. and wa1/ka gt 7.) then stop
   endwhile
   wa_rooti = round(wa_rooti*1.e4)/1.e4
   wrlocs = where(wa_rooti[uniq(wa_rooti)] ne 0. and wa_rooti[uniq(wa_rooti)] lt vae*ka,nwrlocs)

   if wrlocs[0] ne -1 then begin
      mrk = 1
      print,wa_rooti((uniq(wa_rooti))[wrlocs])/ka
      if j eq 0 then begin
         waka_root = wa_rooti((uniq(wa_rooti))[wrlocs])/ka 
         ka_root = replicate(ka,nwrlocs)
      endif else begin
         waka_root = [waka_root,wa_rooti((uniq(wa_rooti))[wrlocs])/ka]
         ka_root = [ka_root,replicate(ka,nwrlocs)]
      endelse
   endif
endfor 

plot,ka_root,waka_root,/xs,/ys,psym=3,yr=[va-va/5,vae+vae/5.]

end
