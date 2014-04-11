
pro init,ro=ro,re=re,vao=vao,vae=vae,co=co,ce=ce,bo=bo,be=be,waka_ini_f=waka_ini_f,waka_ini_ar=waka_ini_ar,nmode=nmode


; ------------------------------------------------------------------
;if ~keyword_set(ro) then begin  
;    print,'init,ro=ro,re=re,vao=vao,vae=vae,co=co,ce=ce,bo=bo,be=be,waka_ini_f=waka_ini_f,waka_ini_ar=waka_ini_ar,nmode=nmode'
;    return
;endif

;--------------------------------------------------------------------

; INPUT:

; nmode = azimuthal wavenumber (n=0: sausage mode, n=1: kink mode)

; OUTPUT:
; ro = internal density
; re = external density
; vao = internal Alfven velocity
; vae = external Alfven velocity
; co = internal sound speed
; ce = external sound speed
; bo = internal magnetic field
; be = external magnetic field
; waka_ini_f = initial guess of w/k of fundamental mode for solving dispersion relation
; waka_ini_ar = initial guesses of w/k of other branches


   dim = 10000
   ka = 4.
   norm = 1.e5
   L = fltarr(dim+1)
   moa2 = fltarr(dim+1)
   mea2 = fltarr(dim+1)
   bslij = fltarr(dim+1)
   bslky = fltarr(dim+1)
   dbslij = fltarr(dim+1)
   dbslky = fltarr(dim+1)
   fracij = fltarr(dim+1)
   fracky = fltarr(dim+1)
   sigi = intarr(dim+1)
   sigk = intarr(dim+1)
   if keyword_set(nmode) eq 1 then n = nmode else n = 0
   mup = 1.25663706*1.e-6
   kboltz = 1.380658*10^(-23.)
   proton = 1.67262158*10^(-27.)
   gamma = 5./3.

; INPUT:
   ne_out_cgs = 3.e9
   ne_in_cgs = 1.e10
   te_out = 2.e6
   te_in = 1.e7
   
   ne_out = ne_out_cgs * 1.e6
   ne_in = ne_in_cgs * 1.e6
   rho_out = ne_out * proton
   rho_in = ne_in * proton
   
; (Uniform magnetic field)
   Bx_out = 0.
   Bz_out = 0.00561431
   By_out = 0.
   Bx_in = 0.
   Bz_in = 0.005
   By_in = 0.
   
   beta = 2*ne_in*te_in*kboltz*mup/(Bx_in^2+By_in^2+Bz_in^2)
   
   ro = rho_in*norm^2
   re = rho_out*norm^2
   vae = sqrt(Bx_out^2+By_out^2+Bz_out^2)/sqrt(mup*re)
   vao = sqrt(Bx_in^2+By_in^2+Bz_in^2)/sqrt(mup*ro)
   co = sqrt(2*gamma*kboltz/proton*te_in)/norm
   ce = sqrt(2*gamma*kboltz/proton*te_out)/norm
   bo = vao*sqrt(mup*ro)
   be = vae*sqrt(mup*re)
   ct = co*vao/sqrt(co^2+vao^2)
   cte = ce*vae/sqrt(ce^2+vae^2)

   vmin = 4.
   va_mx = ceil(max([vao,vae]))
   wa = findgen(dim+1)/float(dim)*(va_mx-vmin)*ka+vmin*ka

   nummoa2 = (ka^2*co^2-wa^2)*(ka^2*vao^2-wa^2)
   dnummoa2 = (ka^2*ct^2-wa^2)*(co^2+vao^2)
   moa2 = nummoa2/dnummoa2
   nummea2 = (ka^2*ce^2-wa^2)*(ka^2*vae^2-wa^2)
   dnummea2 = (ka^2*cte^2-wa^2)*(ce^2+vae^2)
   mea2 = nummea2/dnummea2

   reg1 = where(moa2 gt 0. and mea2 gt 0.)
   reg2 = where(moa2 gt 0. and mea2 lt 0.)
   reg3 = where(moa2 lt 0. and mea2 gt 0.)
   reg4 = where(moa2 lt 0. and mea2 lt 0.)

   if reg1[0] ne -1 then begin
      bslij[reg1] = beseli(sqrt(moa2[reg1]),n,/double)
      bslky[reg1] = beselk(sqrt(mea2[reg1]),n,/double)
      dbslij[reg1] = dbeseli(sqrt(moa2[reg1]),n)
      dbslky[reg1] = dbeselk(sqrt(mea2[reg1]),n)
      fracij[reg1] = dbslij[reg1]/bslij[reg1]
      fracky[reg1] = dbslky[reg1]/bslky[reg1]
      sigi[reg1] = 1.
      sigk[reg1] = 1.
   endif
   if reg2[0] ne -1 then begin
      bslij[reg2] = beseli(sqrt(moa2[reg2]),n,/double)
      bslky[reg2] = besely(sqrt(-mea2[reg2]),n,/double)
      dbslij[reg2] = dbeseli(sqrt(moa2[reg2]),n)
      dbslky[reg2] = dbesely(sqrt(-mea2[reg2]),n)
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
      bslij[reg4] = beselj(sqrt(-moa2[reg4]),n,/double)
      bslky[reg4] = beselk(sqrt(-mea2[reg4]),n,/double)
      dbslij[reg4] = dbeselj(sqrt(-moa2[reg4]),n)
      dbslky[reg4] = dbeselk(sqrt(-mea2[reg4]),n)
      fracij[reg4] = dbslij[reg4]/bslij[reg4]
      fracky[reg4] = dbslky[reg4]/bslky[reg4]
      sigi[reg4] = -1.
      sigk[reg4] = -1.
   endif

   if reg1[0] ne -1 then L[reg1] = wa[reg1]^2*(re*sqrt(sigi[reg1]*moa2[reg1])*dbslij[reg1]*bslky[reg1]-ro*sqrt(sigk[reg1]*mea2[reg1])*dbslky[reg1]*bslij[reg1])-$
    ka^2*(re*vae^2*sqrt(sigi[reg1]*moa2[reg1])*dbslij[reg1]*bslky[reg1]-ro*vao^2*sqrt(sigk[reg1]*mea2[reg1])*dbslky[reg1]*bslij[reg1])
   if reg2[0] ne -1 then L[reg2] = wa[reg2]^2*(re*sqrt(sigi[reg2]*moa2[reg2])*dbslij[reg2]*bslky[reg2]-ro*sqrt(sigk[reg2]*mea2[reg2])*dbslky[reg2]*bslij[reg2])-$
    ka^2*(re*vae^2*sqrt(sigi[reg2]*moa2[reg2])*dbslij[reg2]*bslky[reg2]-ro*vao^2*sqrt(sigk[reg2]*mea2[reg2])*dbslky[reg2]*bslij[reg2])
   if reg3[0] ne -1 then L[reg3] = wa[reg3]^2*(re*sqrt(sigi[reg3]*moa2[reg3])*dbslij[reg3]*bslky[reg3]-ro*sqrt(sigk[reg3]*mea2[reg3])*dbslky[reg3]*bslij[reg3])-$
    ka^2*(re*vae^2*sqrt(sigi[reg3]*moa2[reg3])*dbslij[reg3]*bslky[reg3]-ro*vao^2*sqrt(sigk[reg3]*mea2[reg3])*dbslky[reg3]*bslij[reg3])
   if reg4[0] ne -1 then L[reg4] = wa[reg4]^2*(re*sqrt(sigi[reg4]*moa2[reg4])*dbslij[reg4]*bslky[reg4]-ro*sqrt(sigk[reg4]*mea2[reg4])*dbslky[reg4]*bslij[reg4])-$
    ka^2*(re*vae^2*sqrt(sigi[reg4]*moa2[reg4])*dbslij[reg4]*bslky[reg4]-ro*vao^2*sqrt(sigk[reg4]*mea2[reg4])*dbslky[reg4]*bslij[reg4])

bound = 5e-2

   window,0
;plot,wa/ka,L/max(sqrt(abs(dnummoa2*dnummea2))),/xs,/ys,psym=3;,xr=[2,5],yr=[-1,1]
;locs = where(abs(L/max(sqrt(abs(dnummoa2*dnummea2)))) lt 1.e-5)

   plot,wa/ka,L,/xs,/ys,psym=3,yr=[-15,15]  ;,yr=[-0.01,0.01]
   locs = where(abs(L) lt bound)
   if locs[0] ne -1 then begin
      ;waka_ini_f = max(wa[locs]/ka[locs])             ; Two possibilities for setting an initial guess
      ;waka_ini_f = wa(where(L eq min(abs(L))))/ka     ; Normally each one selects a different branch to compute roots
      waka_ini_f = 14.3787
      waka_ini_0 = round(wa[locs]/ka[locs]*1.e5)/1.e5
      wrlocs = where(waka_ini_0[uniq(waka_ini_0)] gt vao*1.01)
      waka_ini_ar = waka_ini_f                                     ; Copied from GS pressure balance file (normally via waka_ini_f and wrlocs)
      ;waka_ini_ar = waka_ini_0[(uniq(waka_ini_0))[wrlocs]]
      print,'waka_ini_f = ', waka_ini_f
   endif else begin
      print,'no solution'
   endelse

end
   
