
pro vel_modes, waka_root=waka_root, ka_root=ka_root, dimx=dimx, reg0=reg0, reg1=reg1, reg2=reg2, reg3=reg3, reg4=reg4, vr_md=vr_md, vt_md=vt_md, vz_md=vz_md, pr_md=pr_md, ptot_md=ptot_md, rr_md=rr_md, br_md=br_md, bt_md=bt_md, bz_md=bz_md, btot_md=btot_md, gridr=gridr, gridx=gridx, aa=aa,nmode=nmode

if ~keyword_set(waka_root) then begin
   print,'vel_modes, waka_root=waka_root, ka_root=ka_root, dimx=dimx, reg0=reg0, reg1=reg1, reg2=reg2, reg3=reg3, reg4=reg4, vr_md=vr_md, vt_md=vt_md, vz_md=vz_md, pr_md=pr_md, ptot_md=ptot_md, rr_md=rr_md, br_md=br_md, bt_md=bt_md, bz_md=bz_md, btot_md=btot_md, gridr=gridr, gridx=gridx, aa=aa,nmode=nmode'
   return
endif

; Calculates the initial setup of the modulation of the MHD mode
; on thermodynamic and geometrical quantities (vr, vt, vz, pr, rr) 
; according to Edwin & Roberts 1983
; INPUT:
; waka_root = w/k roots of dispersion relation
; ka_root = values of ka corresponding to w/k roots (a = radius of cylinder)
; dimx = dimension along x axis
; OUTPUT:
; regi (i=0,...,4) = positions in ka_root array corresponding to
; regions of different behavior (trapped, leaky,...)
; gridr = grid along r axis (radial direction)
; vr_md = vr(k,x) : array of radial velocity (wavelength, x position)
; vt_md = vt(k,x) : array of azimuthal velocity (wavelength, x position)
; vz_md = vz(k,x) : array of longitudinal velocity (wavelength, x position)
; pr_md = pr(k,x) : array of pressure perturbation along radial direction (x axis)
; rr_md = rr(k,x) : array of density perturbation along radial direction (x axis)
; br_md = rr(k,x) : array of radial magnetic field perturbation along radial direction (x axis)
; bt_md = rr(k,x) : array of azimuthal magnetic field perturbation along radial direction (x axis)
; bz_md = rr(k,x) : array of longitudinal magnetic field perturbation along radial direction (x axis)

common vars1, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, n, ka

if keyword_set(nmode) then n = nmode else n = 0 ; default is sausage mode

co = A1 & vao = A2 & ct = A3 & ro = A4 & bo = A5
ce = B1 & vae = B2 & cte = B3 & re = B4 & be = B5

aa = 10.
gridx = findgen(dimx)/float(dimx-1)*aa*4.
;if n eq 0 then gridx = findgen(dimx)/float(dimx-1)*50. else gridx = findgen(dimx)/float(dimx-1)*100.
r0 = gridx[n_elements(gridx)-1]/2.
;if n eq 0 then begin
;   gridr = gridx-r0
;   dimr = dimx
;endif else begin
;dimr = dimx
;hypo = sqrt(2)*abs(r0-gridx[dimx-1])
resx = gridx[dimx-1]/dimx
;dimr = hypo/resx*2
dimr = dimx
gridr = findgen(dimr)/(dimr-1)*abs(r0-gridx[dimx-1])
;gridrp = findgen(dimr)/(dimr-1)*hypo*2
;gridr = gridrp-gridrp[dimr-1]/2.
;   gridr = gridx[dimx/2:dimx-1]-r0
;endelse

polyind = 5./3.
dim = n_elements(waka_root)
bslij = fltarr(dim,dimr)
bslky = fltarr(dim,dimr)
vr_md = fltarr(dim,dimr)
vt_md = fltarr(dim,dimr)
vz_md = fltarr(dim,dimr)
pr_md = fltarr(dim,dimr)
ptot_md = fltarr(dim,dimr)
rr_md = fltarr(dim,dimr)
br_md = fltarr(dim,dimr)
bt_md = fltarr(dim,dimr)
bz_md = fltarr(dim,dimr)
btot_md = fltarr(dim,dimr)
mor2 = fltarr(dim,dimr)
mer2 = fltarr(dim,dimr)
R = fltarr(dim,dimr)
sigi = fltarr(dim)
sigk = fltarr(dim)
aa1 = fltarr(dim)
aa0 = fltarr(dim)
aa0[*] = 1.

nummor2 = (co^2-waka_root^2)*(vao^2-waka_root^2)
dnummor2 = (ct^2-waka_root^2)*(co^2+vao^2)
nummer2 = (ce^2-waka_root^2)*(vae^2-waka_root^2)
dnummer2 = (cte^2-waka_root^2)*(ce^2+vae^2)
moa2 = ka_root^2*nummor2/dnummor2
mea2 = ka_root^2*nummer2/dnummer2

locout = where(abs(gridr) ge aa,nlocout)
locin = where(abs(gridr) lt aa,nlocin)

;neglocin = where(gridr[locin] lt 0.)
;poslocin = where(gridr[locin] ge 0.)
;neglocout = where(gridr[locout] lt 0.)
;poslocout = where(gridr[locout] ge 0.)

;signr = intarr(dimr)
;signr[locin[poslocin]] = 1
;signr[locin[neglocin]] = -1
;signr[locout[poslocout]] = 1
;signr[locout[neglocout]] = -1

thrs = 1.e-4
reg0 = where(abs(moa2) lt thrs or abs(mea2) lt thrs)
reg1 = where(moa2 gt thrs and mea2 gt thrs)
reg2 = where(moa2 gt thrs and mea2 lt -thrs)
reg3 = where(moa2 lt -thrs and mea2 gt thrs)
reg4 = where(moa2 lt -thrs and mea2 lt -thrs)

if reg1[0] ne -1 then begin
   sigi[reg1] = 1.
   sigk[reg1] = 1.
   aa1[reg1] = aa0[reg1]*beseli(sqrt(sigi[reg1]*moa2[reg1]),n,/double)/beselk(sqrt(sigk[reg1]*mea2[reg1]),n,/double)*(co^2+vao^2)*(waka_root[reg1]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg1]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg1,locin[i]] = moa2[reg1]*(gridr[locin[i]]/aa)^2
      vr_md[reg1,locin[i]] = (vao^2+co^2)*(waka_root[reg1]^2-ct^2)/(waka_root[reg1]^2*(waka_root[reg1]^2-vao^2))*aa0[reg1]*sqrt(sigi[reg1]*moa2[reg1]/aa^2)*dbeseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n)*(aa/ka_root[reg1])^2
      if gridr[locin[i]] eq 0 then begin         
         vt_md[reg1,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg1]^2)*aa0[reg1]*n*1./((waka_root[reg1]*ka_root[reg1])^2-vao^2*ka_root[reg1]^2)*aa^2*(1./gamma(n+1)*(sqrt(sigi[reg1]*moa2[reg1]/aa)*0.5)^n)
      endif else begin
         vt_md[reg1,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg1]^2)*aa0[reg1]*n*beseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n,/double)/((waka_root[reg1]*ka_root[reg1])^2-vao^2*ka_root[reg1]^2)/gridr[locin[i]]*aa^2
      endelse
      vz_md[reg1,locin[i]] = 1./waka_root[reg1]^2/ka_root[reg1]*aa*co^2*aa0[reg1]*beseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n,/double)
      pr_md[reg1,locin[i]] = co^2*ro/waka_root[reg1]/ka_root[reg1]*aa*aa0[reg1]*beseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n,/double)
      rr_md[reg1,locin[i]] = ro/waka_root[reg1]/ka_root[reg1]*aa*aa0[reg1]*beseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n,/double)
      br_md[reg1,locin[i]] = bo/waka_root[reg1]*vr_md[reg1,locin[i]]
      bt_md[reg1,locin[i]] = bo/waka_root[reg1]*vt_md[reg1,locin[i]]
      bz_md[reg1,locin[i]] = bo/waka_root[reg1]/ka_root[reg1]*aa*(1.-co^2/waka_root[reg1]^2)*aa0[reg1]*beseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n,/double)
      btot_md[reg1,locin[i]] = sqrt(br_md[reg1,locin[i]]^2+bt_md[reg1,locin[i]]^2+bz_md[reg1,locin[i]]^2+2*bo*bz_md[reg1,locin[i]]+bo^2)
      ptot_md[reg1,locin[i]] = (co^2+vao^2-co^2*vao^2/waka_root[reg1]^2)*ro/waka_root[reg1]/ka_root[reg1]*aa*aa0[reg1]*beseli(sqrt(sigi[reg1]*mor2[reg1,locin[i]]),n,/double)
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg1,locout[i]] = mea2[reg1]*(gridr[locout[i]]/aa)^2
      vr_md[reg1,locout[i]] = (vae^2+ce^2)*(waka_root[reg1]^2-cte^2)/(waka_root[reg1]^2*(waka_root[reg1]^2-vae^2))*aa1[reg1]*sqrt(sigk[reg1]*mea2[reg1]/aa^2)*dbeselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n)*(aa/ka_root[reg1])^2
      vt_md[reg1,locout[i]] = -(vae^2+ce^2-vae^2*ce^2/waka_root[reg1]^2)*aa1[reg1]*n*beselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n,/double)/((waka_root[reg1]*ka_root[reg1])^2-vae^2*ka_root[reg1]^2)/gridr[locout[i]]*aa^2
      vz_md[reg1,locout[i]] = 1./waka_root[reg1]^2/ka_root[reg1]*aa*ce^2*aa1[reg1]*beselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n,/double)
      pr_md[reg1,locout[i]] = ce^2*re/waka_root[reg1]/ka_root[reg1]*aa*aa1[reg1]*beselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n,/double)
      rr_md[reg1,locout[i]] = re/waka_root[reg1]/ka_root[reg1]*aa*aa1[reg1]*beselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n,/double)
      br_md[reg1,locout[i]] = be/waka_root[reg1]*vr_md[reg1,locout[i]]
      bt_md[reg1,locout[i]] = be/waka_root[reg1]*vt_md[reg1,locout[i]]
      bz_md[reg1,locout[i]] = be/waka_root[reg1]/ka_root[reg1]*aa*(1.-ce^2/waka_root[reg1]^2)*aa1[reg1]*beselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n,/double)
      btot_md[reg1,locout[i]] = sqrt(br_md[reg1,locout[i]]^2+bt_md[reg1,locout[i]]^2+bz_md[reg1,locout[i]]^2+2*be*bz_md[reg1,locout[i]]+be^2)
      ptot_md[reg1,locout[i]] = (ce^2+vae^2-ce^2*vae^2/waka_root[reg1]^2)*re/waka_root[reg1]/ka_root[reg1]*aa*aa1[reg1]*beselk(sqrt(sigk[reg1]*mer2[reg1,locout[i]]),n,/double)
   endfor
endif

if reg2[0] ne -1 then begin
   sigi[reg2] = 1.
   sigk[reg2] = -1.
   aa1[reg2] = aa0[reg2]*beseli(sqrt(sigi[reg2]*moa2[reg2]),n,/double)/besely(sqrt(sigk[reg2]*mea2[reg2]),n,/double)*(co^2+vao^2)*(waka_root[reg2]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg2]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg2,locin[i]] = moa2[reg2]*(gridr[locin[i]]/aa)^2
      vr_md[reg2,locin[i]] = (vao^2+co^2)*(waka_root[reg2]^2-ct^2)/(waka_root[reg2]^2*(waka_root[reg2]^2-vao^2))*aa0[reg2]*sqrt(sigi[reg2]*moa2[reg2]/aa^2)*dbeseli(sqrt(sigi[reg2]*mor2[reg2,locin[i]]),n)*(aa/ka_root[reg2])^2
      if gridr[locin[i]] eq 0 then begin         
         vt_md[reg2,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg2]^2)*aa0[reg2]*n*1./((waka_root[reg2]*ka_root[reg2])^2-vao^2*ka_root[reg2]^2)*aa^2*(1./gamma(n+1)*(sqrt(sigi[reg2]*moa2[reg2]/aa)*0.5)^n)
      endif else begin
         vt_md[reg2,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg2]^2)*aa0[reg2]*n*beseli(sqrt(sigi[reg2]*mor2[reg2,locin[i]]),n,/double)/((waka_root[reg2]*ka_root[reg2])^2-vao^2*ka_root[reg2]^2)/gridr[locin[i]]*aa^2
      endelse
      vz_md[reg2,locin[i]] = 1./waka_root[reg2]^2/ka_root[reg2]*aa*co^2*aa0[reg2]*beseli(sqrt(sigi[reg2]*mor2[reg2,locin[i]]),n,/double)
      pr_md[reg2,locin[i]] = co^2*ro/waka_root[reg2]/ka_root[reg2]*aa*aa0[reg2]*beseli(sigi[reg2]*sqrt(mor2[reg2,locin[i]]),n,/double)
      rr_md[reg2,locin[i]] = ro/waka_root[reg2]/ka_root[reg2]*aa*aa0[reg2]*beseli(sqrt(sigi[reg2]*mor2[reg2,locin[i]]),n,/double)
      br_md[reg2,locin[i]] = bo/waka_root[reg2]*vr_md[reg2,locin[i]]
      bt_md[reg2,locin[i]] = bo/waka_root[reg2]*vt_md[reg2,locin[i]]
      bz_md[reg2,locin[i]] = bo/waka_root[reg2]/ka_root[reg2]*aa*(1.-co^2/waka_root[reg2]^2)*aa0[reg2]*beseli(sqrt(sigi[reg2]*mor2[reg2,locin[i]]),n,/double)
      btot_md[reg2,locin[i]] = sqrt(br_md[reg2,locin[i]]^2+bt_md[reg2,locin[i]]^2+bz_md[reg2,locin[i]]^2+2*bo*bz_md[reg2,locin[i]]+bo^2)
      ptot_md[reg2,locin[i]] = (co^2+vao^2-co^2*vao^2/waka_root[reg2]^2)*ro/waka_root[reg2]/ka_root[reg2]*aa*aa0[reg2]*beseli(sqrt(sigi[reg2]*mor2[reg2,locin[i]]),n,/double)
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg2,locout[i]] = mea2[reg2]*(gridr[locout[i]]/aa)^2
      vr_md[reg2,locout[i]] = (vae^2+ce^2)*(waka_root[reg2]^2-cte^2)/(waka_root[reg2]^2*(waka_root[reg2]^2-vae^2))*aa1[reg2]*sqrt(sigk[reg2]*mea2[reg2]/aa^2)*dbesely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n)*(aa/ka_root[reg2])^2
      vt_md[reg2,locout[i]] = -(vae^2+ce^2-vae^2*ce^2/waka_root[reg2]^2)*aa1[reg2]*n*besely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n,/double)/((waka_root[reg2]*ka_root[reg2])^2-vae^2*ka_root[reg2]^2)/gridr[locout[i]]*aa^2
      vz_md[reg2,locout[i]] = 1./waka_root[reg2]^2/ka_root[reg2]*aa*ce^2*aa1[reg2]*besely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n,/double)
      pr_md[reg2,locout[i]] = ce^2*re/waka_root[reg2]/ka_root[reg2]*aa*aa1[reg2]*besely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n,/double)
      rr_md[reg2,locout[i]] = re/waka_root[reg2]/ka_root[reg2]*aa*aa1[reg2]*besely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n,/double)
      br_md[reg2,locout[i]] = be/waka_root[reg2]*vr_md[reg2,locout[i]]
      bt_md[reg2,locout[i]] = be/waka_root[reg2]*vt_md[reg2,locout[i]]
      bz_md[reg2,locout[i]] = be/waka_root[reg2]/ka_root[reg2]*aa*(1.-ce^2/waka_root[reg2]^2)*aa1[reg2]*besely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n,/double)
      btot_md[reg2,locout[i]] = sqrt(br_md[reg2,locout[i]]^2+bt_md[reg2,locout[i]]^2+bz_md[reg2,locout[i]]^2+2*be*bz_md[reg2,locout[i]]+be^2)
      ptot_md[reg2,locout[i]] = (ce^2+vae^2-ce^2*vae^2/waka_root[reg2]^2)*re/waka_root[reg2]/ka_root[reg2]*aa*aa1[reg2]*besely(sqrt(sigk[reg2]*mer2[reg2,locout[i]]),n,/double)
   endfor
endif

if reg3[0] ne -1 then begin
   sigi[reg3] = -1.
   sigk[reg3] = 1.
;   aa1[reg3] = aa0[reg3]*beselj(sqrt(sigi[reg3]*moa2[reg3]/aa^2*gridr[locin[-1]]^2),n,/double)/beselk(sqrt(sigk[reg3]*mea2[reg3]/aa^2*gridr[locout[0]]^2),n,/double)*(co^2+vao^2)*(waka_root[reg3]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg3]^2-cte^2))*ro/re
   aa1[reg3] = aa0[reg3]*beselj(sqrt(sigi[reg3]*moa2[reg3]),n,/double)/beselk(sqrt(sigk[reg3]*mea2[reg3]),n,/double)*(co^2+vao^2)*(waka_root[reg3]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg3]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg3,locin[i]] = moa2[reg3]*(gridr[locin[i]]/aa)^2
      vr_md[reg3,locin[i]] = (vao^2+co^2)*(waka_root[reg3]^2-ct^2)/(waka_root[reg3]^2*(waka_root[reg3]^2-vao^2))*aa0[reg3]*sqrt(sigi[reg3]*moa2[reg3]/aa^2)*dbeselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n)*(aa/ka_root[reg3])^2
      if gridr[locin[i]] eq 0 then begin         
         if n gt 0 then vt_md[reg3,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg3]^2)*aa0[reg3]*n*1./((waka_root[reg3]*ka_root[reg3])^2-vao^2*ka_root[reg3]^2)*aa^2*1./gamma(n+1)*(sqrt(sigi[reg3]*moa2[reg3]/aa^2)*0.5)^n*(gridr[locin[i]])^(n-1) else vt_md[reg3,locin[i]] = 0.
      endif else begin
         vt_md[reg3,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg3]^2)*aa0[reg3]*n*beselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n,/double)/((waka_root[reg3]*ka_root[reg3])^2-vao^2*ka_root[reg3]^2)/gridr[locin[i]]*aa^2
      endelse
      vz_md[reg3,locin[i]] = 1./waka_root[reg3]^2/ka_root[reg3]*aa*co^2*aa0[reg3]*beselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n,/double)
      pr_md[reg3,locin[i]] = co^2*ro/waka_root[reg3]/ka_root[reg3]*aa*aa0[reg3]*beselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n,/double)
      rr_md[reg3,locin[i]] = ro/waka_root[reg3]/ka_root[reg3]*aa*aa0[reg3]*beselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n,/double)
      br_md[reg3,locin[i]] = bo/waka_root[reg3]*vr_md[reg3,locin[i]]
      bt_md[reg3,locin[i]] = bo/waka_root[reg3]*vt_md[reg3,locin[i]]
      bz_md[reg3,locin[i]] = bo/waka_root[reg3]/ka_root[reg3]*aa*(1.-co^2/waka_root[reg3]^2)*aa0[reg3]*beselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n,/double)
      btot_md[reg3,locin[i]] = sqrt(br_md[reg3,locin[i]]^2+bt_md[reg3,locin[i]]^2+bz_md[reg3,locin[i]]^2+2*bo*bz_md[reg3,locin[i]]+bo^2)
      ptot_md[reg3,locin[i]] = (co^2+vao^2-co^2*vao^2/waka_root[reg3]^2)*ro/waka_root[reg3]/ka_root[reg3]*aa*aa0[reg3]*beselj(sqrt(sigi[reg3]*mor2[reg3,locin[i]]),n,/double)
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg3,locout[i]] = mea2[reg3]*(gridr[locout[i]]/aa)^2
      vr_md[reg3,locout[i]] = (vae^2+ce^2)*(waka_root[reg3]^2-cte^2)/(waka_root[reg3]^2*(waka_root[reg3]^2-vae^2))*aa1[reg3]*sqrt(sigk[reg3]*mea2[reg3]/aa^2)*dbeselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n)*(aa/ka_root[reg3])^2
      vt_md[reg3,locout[i]] = -(vae^2+ce^2-vae^2*ce^2/waka_root[reg3]^2)*aa1[reg3]*n*beselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n,/double)/((waka_root[reg3]*ka_root[reg3])^2-vae^2*ka_root[reg3]^2)/gridr[locout[i]]*aa^2
      vz_md[reg3,locout[i]] = 1./waka_root[reg3]^2/ka_root[reg3]*aa*ce^2*aa1[reg3]*beselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n,/double)
      pr_md[reg3,locout[i]] = ce^2*re/waka_root[reg3]/ka_root[reg3]*aa*aa1[reg3]*beselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n,/double)
      rr_md[reg3,locout[i]] = re/waka_root[reg3]/ka_root[reg3]*aa*aa1[reg3]*beselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n,/double)
      br_md[reg3,locout[i]] = be/waka_root[reg3]*vr_md[reg3,locout[i]]
      bt_md[reg3,locout[i]] = be/waka_root[reg3]*vt_md[reg3,locout[i]]
      bz_md[reg3,locout[i]] = be/waka_root[reg3]/ka_root[reg3]*aa*(1.-ce^2/waka_root[reg3]^2)*aa1[reg3]*beselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n,/double)
      btot_md[reg3,locout[i]] = sqrt(br_md[reg3,locout[i]]^2+bt_md[reg3,locout[i]]^2+bz_md[reg3,locout[i]]^2+2*be*bz_md[reg3,locout[i]]+be^2)
      ptot_md[reg3,locout[i]] = (ce^2+vae^2-ce^2*vae^2/waka_root[reg3]^2)*re/waka_root[reg3]/ka_root[reg3]*aa*aa1[reg3]*beselk(sqrt(sigk[reg3]*mer2[reg3,locout[i]]),n,/double)
   endfor
endif

if reg4[0] ne -1 then begin
   sigi[reg4] = -1.
   sigk[reg4] = -1.
   aa1[reg4] = aa0[reg4]*beselj(sqrt(sigi[reg4]*moa2[reg4]),n,/double)/besely(sqrt(sigk[reg4]*mea2[reg4]),n,/double)*(co^2+vao^2)*(waka_root[reg4]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg4]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg4,locin[i]] = moa2[reg4]*(gridr[locin[i]]/aa)^2
      vr_md[reg4,locin[i]] = (vao^2+co^2)*(waka_root[reg4]^2-ct^2)/(waka_root[reg4]^2*(waka_root[reg4]^2-vao^2))*aa0[reg4]*sqrt(sigi[reg4]*moa2[reg4]/aa^2)*dbeselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n)*(aa/ka_root[reg4])^2
      if gridr[locin[i]] eq 0 then begin         
         vt_md[reg4,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg4]^2)*aa0[reg4]*n*1./((waka_root[reg4]*ka_root[reg4])^2-vao^2*ka_root[reg4]^2)*aa^2*(1./gamma(n+1)*(sqrt(sigi[reg4]*moa2[reg4]/aa)*0.5)^n)
      endif else begin
         vt_md[reg4,locin[i]] = -(vao^2+co^2-vao^2*co^2/waka_root[reg4]^2)*aa0[reg4]*n*beselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n,/double)/((waka_root[reg4]*ka_root[reg4])^2-vao^2*ka_root[reg4]^2)/gridr[locin[i]]*aa^2
      endelse
      vz_md[reg4,locin[i]] = 1./waka_root[reg4]^2/ka_root[reg4]*aa*co^2*aa0[reg4]*beselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n,/double)
      pr_md[reg4,locin[i]] = co^2*ro/waka_root[reg4]/ka_root[reg4]*aa*aa0[reg4]*beselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n,/double)
      rr_md[reg4,locin[i]] = ro/waka_root[reg4]/ka_root[reg4]*aa*aa0[reg4]*beselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n,/double)
      br_md[reg4,locin[i]] = bo/waka_root[reg4]*vr_md[reg4,locin[i]]
      bt_md[reg4,locin[i]] = bo/waka_root[reg4]*vt_md[reg4,locin[i]]
      bz_md[reg4,locin[i]] = bo/waka_root[reg4]/ka_root[reg4]*aa*(1.-co^2/waka_root[reg4]^2)*aa0[reg4]*beselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n,/double)
      btot_md[reg4,locin[i]] = sqrt(br_md[reg4,locin[i]]^2+bt_md[reg4,locin[i]]^2+bz_md[reg4,locin[i]]^2+2*bo*bz_md[reg4,locin[i]]+bo^2)
      ptot_md[reg4,locin[i]] = (co^2+vao^2-co^2*vao^2/waka_root[reg4]^2)*ro/waka_root[reg4]/ka_root[reg4]*aa*aa0[reg4]*beselj(sqrt(sigi[reg4]*mor2[reg4,locin[i]]),n,/double)
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg4,locout[i]] = mea2[reg4]*(gridr[locout[i]]/aa)^2
      vr_md[reg4,locout[i]] = (vae^2+ce^2)*(waka_root[reg4]^2-cte^2)/(waka_root[reg4]^2*(waka_root[reg4]^2-vae^2))*aa1[reg4]*sqrt(sigk[reg4]*mea2[reg4]/aa^2)*dbesely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n)*(aa/ka_root[reg4])^2
      vt_md[reg4,locout[i]] = -(vae^2+ce^2-vae^2*ce^2/waka_root[reg4]^2)*aa1[reg4]*n*besely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n,/double)/((waka_root[reg4]*ka_root[reg4])^2-vae^2*ka_root[reg4]^2)/gridr[locout[i]]*aa^2
      vz_md[reg4,locout[i]] = 1./waka_root[reg4]^2/ka_root[reg4]*aa*ce^2*aa1[reg4]*besely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n,/double)
      pr_md[reg4,locout[i]] = ce^2*re/waka_root[reg4]/ka_root[reg4]*aa*aa1[reg4]*besely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n,/double)
      rr_md[reg4,locout[i]] = re/waka_root[reg4]/ka_root[reg4]*aa*aa1[reg4]*besely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n,/double)
      br_md[reg4,locout[i]] = be/waka_root[reg4]*vr_md[reg4,locout[i]]
      bt_md[reg4,locout[i]] = be/waka_root[reg4]*vt_md[reg4,locout[i]]
      bz_md[reg4,locout[i]] = be/waka_root[reg4]/ka_root[reg4]*aa*(1.-ce^2/waka_root[reg4]^2)*aa1[reg4]*besely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n,/double)
      btot_md[reg4,locout[i]] = sqrt(br_md[reg4,locout[i]]^2+bt_md[reg4,locout[i]]^2+bz_md[reg4,locout[i]]^2+2*bo*bz_md[reg4,locout[i]]+bo^2)
      ptot_md[reg4,locout[i]] = (ce^2+vae^2-ce^2*vae^2/waka_root[reg4]^2)*re/waka_root[reg4]/ka_root[reg4]*aa*aa1[reg4]*besely(sqrt(sigk[reg4]*mer2[reg4,locout[i]]),n,/double)
   endfor
endif

;print, size(vz_md)

end
