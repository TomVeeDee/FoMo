
PRO tubemodes, ro=ro, re=re, vao=vao, vae=vae, co=co, ce=ce, bo=bo, be=be, waka_ini_f=waka_ini_f, waka_root=waka_root, ka_root=ka_root, nmode = nmode, normcase = normcase

if ~keyword_set(ro) then begin
   print,'tubemodes, ro=ro, re=re, vao=vao, vae=vae, co=co, ce=ce, bo=bo, be=be, waka_ini_f=waka_ini_f, waka_root=waka_root, ka_root=ka_root, nmode = nmode, normcase = normcase'
   return
endif

; Solves dispersion relation based on parameters ro, re, vao, vae, co,
; ce, waka_ini_f.
; INPUT:
; ro = internal density: rho[dimx/2,dimy/2,dimz/2]
; re = external density: rho[0,0,0]
; vao = internal Alfven velocity: calculated at (dimx/2,dimy/2,dimz/2)
; vae = external Alfven velocity: calculated at (0,0,0)
; co = internal sound speed: calculated at (dimx/2,dimy/2,dimz/2)
; ce = external sound speed: calculated at (0,0,0)
; bo = internal magnetic field calculated at (dimx/2,dimy/2,dimz/2)
; be = external magnetic field calculated at (0,0,0)
; waka_ini_f = initial guess of w/k for solving dispersion relation
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

gamma = 5./3.
nka = 400 ; number of roots
ka_ar = findgen(nka)/(nka-1.)*3.99+0.01 ; initial values for ka
ka_root = fltarr(nka)
waka_root = fltarr(nka)
if keyword_set(normcase) then begin
   co = 1.
   ce = 0.5 * co
   vao = 2 * co
   vae = 5 * co
   ro = 1.
   re = 0.1
   bo = vao*sqrt(ro)
   be = vae*sqrt(re)
endif

ct = co*vao/sqrt(co^2+vao^2)
cte = ce*vae/sqrt(ce^2+vae^2)

A1 = co & A2 = vao & A3 = ct & A4 = ro & A5 = bo
B1 = ce & B2 = vae & B3 = cte & B4 = re & B5 = be
mrk = 0 & mrk2 = 0
waka_j = waka_ini_f

for j=0,nka-1 do begin
   ka_arrev = reverse(ka_ar)
   ka = ka_arrev[j]
   wa0 = waka_j*ka
   wa_rooti = 0.
   dwa1 = 1.
   k = 0
   if max(waka_root) ne 0 then begin 
      wkmx = ([max(waka_root),!c])[0]
      wkmxlc = ([max(waka_root),!c])[1]
      wa0 = wkmx*ka_root[wkmxlc] 
   endif
   print,'inital=',wa0
   while(dwa1 gt 1.e-1) do begin
      stepg = 0.001
      wa = wa0*(1.-stepg*k)
      wa_2 = (wa0*(1.-stepg*k)-0.05*ka)>(vao*1.01*ka)
      wa_3 = wa0*(1.-stepg*k)+0.05*ka
      iX1 = [wa,wa_2,wa_3]
      wa1 = FX_ROOT(iX1,'wkeqsyst',/double,tol=1.e-4)
      if k eq 0 then begin
         wa_rooti = wa1 
         k = k+1
      endif else begin
         dwa1 = wa/ka-vao
         k = k+1
         if wa1/ka gt 0.98*vae then begin
            waka_j=vao*1.1
            mrk = 1
         endif         
         if imaginary(wa1) eq 0. and wa1 gt 0. then begin 
            wa_rooti = [wa_rooti,wa1]
         endif
      endelse
   endwhile
   wa_rooti = round(wa_rooti*1.e4)/1.e4
   wrlocs = where(wa_rooti[uniq(wa_rooti)] ne 0. and wa_rooti[uniq(wa_rooti)] lt vae*ka,nwrlocs)
   
   if wrlocs[0] ne -1 then begin
      mrk = 1 & mrk2 = 1
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

plot,ka_root,waka_root,/xs,/ys,psym=3,yr=[vao-vao/5,vae+vae/5.],xr=[0,4]

end
