
PRO tubemodes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, bo_int=bo, be_ext=be, waka_ini=waka_ini, waka_root, ka_root, nmode = nmode, normcase = normcase

if n_params(0) lt 1 then begin
   print,'tubemodes, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, bo_int=bo, be_ext=be, waka_ini=waka_ini, waka_root, ka_root [,normcase = normcase]'
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
; bo = internal magnetic field calculated at (dimx/2,dimy/2,dimz/2)
; be = external magnetic field calculated at (0,0,0)
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
common vars1, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, n, ka

if keyword_set(nmode) then n = nmode else n = 0 ; default is sausage mode

;xm = round(dimx/2.)
;ym = round(dimy/2.)
;zm = round(dimz/2.)

gamma = 5./3.
mup=1.25663706*1.e-6
nka = 1000 ; number of roots
ka_ar = findgen(nka)/(nka-1.)*4. ; initial values for ka
ka_root = fltarr(nka)
waka_root = fltarr(nka)

if keyword_set(normcase) then begin
   co = 1.
   ce = 0.5 * co
   va = 2 * co
   vae = 5 * co
   ro = 1.
   re = 0.1
   bo = va*sqrt(ro)
   be = vae*sqrt(re)
   ck = sqrt((ro*va^2+re*vae^2)/(ro+re))
endif

ct = co*va/sqrt(co^2+va^2)
cte = ce*vae/sqrt(ce^2+vae^2)

A1 = co & A2 = va & A3 = ct & A4 = ro & A5 = bo
B1 = ce & B2 = vae & B3 = cte & B4 = re & B5 = be
mrk = 0 & mrk2 = 0
waka_j = waka_ini;
dwa1 = 1.
j = 0.
;for j=1,nka-1 do begin
 while(dwa1 gt 1.e-1) do begin
   ka = (reverse(ka_ar))[j]
   waka0 = waka_j
   wk_rooti = 0.
   if max(waka_root) ne 0 then begin
      wkmx = ([max(waka_root),!c])[0]
      wkmxlc = ([max(waka_root),!c])[1]
      waka0 = wkmx
   endif
print,'inital=',waka0

;   while(dwa1 gt 1.e-3) do begin
      stepg = 0.001
      waka = waka0*(1.+stepg)
      waka_2 = waka0*(1.-stepg);+0.1
      waka_3 = waka0*(1.+stepg*3);+0.2
      iX1 = [waka,waka_2,waka_3]*ka
      waka1 = FX_ROOT(iX1,'wkeqsyst',/double,tol=1.e-4)/ka
;      waka = waka0*(1.+stepg)
;      waka_2 = waka0*(1.-stepg*2);+0.1
;      waka_3 = waka0*(1.+stepg*4);+0.2
;      iX2 = [waka,waka_2,waka_3]*ka
;      waka2 = FX_ROOT(iX2,'wkeqsyst',/double,tol=1.e-5)/ka
      if imaginary(waka1) eq 0. then begin
         print,ka,waka1
         wk_rooti = [wk_rooti,waka1]
         dwa1 = vae - waka1
      endif
;   endwhile
   wk_rooti = round(wk_rooti*1.e5)/1.e5
   wrlocs = where(wk_rooti[uniq(wk_rooti)] ne 0. and wk_rooti[uniq(wk_rooti)] lt vae,nwrlocs)

   if wrlocs[0] ne -1 then begin
      mrk = 1 & mrk2 = 1
      print,wk_rooti((uniq(wk_rooti))[wrlocs])
      if j eq 0 then begin
         waka_root = wk_rooti((uniq(wk_rooti))[wrlocs])
         ka_root = replicate(ka,nwrlocs)
      endif else begin
         waka_root = [waka_root,wk_rooti((uniq(wk_rooti))[wrlocs])]
         ka_root = [ka_root,replicate(ka,nwrlocs)]
      endelse
   endif
   j = j+1
endwhile
;endfor 

plot,ka_root,waka_root,/xs,/ys,psym=3,yr=[va-va/5,vae+vae/5.]

end
