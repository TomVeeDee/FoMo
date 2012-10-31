
pro vel_modes, waka_root=waka_root, ka_root=ka_root, dimx=dimx, dimy=dimy, reg0, reg1, reg2, reg3, reg4, gridx, vr_md, vt_md, vz_md, pr_md, rr_md

if n_params(0) lt 1 then begin
   print,'vel_modes, waka_root=waka_root, ka_root=ka_root, dimx=dimx, dimy=dimy, reg0, reg1, reg2, reg3, reg4, gridx, vr_md, vt_md, vz_md, pr_md, rr_md'
   return
endif

; Calculates the initial setup of the modulation of the MHD mode
; on thermodynamic and geometrical quantities (vr, vt, vz, pr, rr) 
; according to Edwin & Roberts 1983
; INPUT:
; waka_root = w/k roots of dispersion relation
; ka_root = values of ka corresponding to w/k roots (a = radius of cylinder)
; dimx = dimension along x axis of structure
; dimy = dimension along y axis of structure
; OUTPUT:
; regi (i=0,...,4) = positions in ka_root array corresponding to
; regions of different behavior (trapped, leaky,...)
; gridx = grid along x axis (perpendicular to cylinder axis)
; vr_md = vr(k,x) : array of radial velocity (wavelength, x position)
; vt_md = vt(k,x) : array of azimuthal velocity (wavelength, x position)
; vz_md = vz(k,x) : array of longitudinal velocity (wavelength, x position)
; pr_md = pr(k,x) : array of pressure perturbation along radial direction (x axis)
; rr_md = rr(k,x) : array of density perturbation along radial direction (x axis)

common vars1, A1, A2, A3, A4, B1, B2, B3, B4, n, ka

co = A1 & va = A2 & ct = A3 & ro = A4
ce = B1 & vae = B2 & cte = B3 & re = B4

gridx = findgen(dimx)/float(dimx-1)*50.
gridy = findgen(dimy)/float(dimy-1)*50.

dimx = n_elements(gridx)
dimy = n_elements(gridy)

n = 0
r0=gridx[-1]/2.
aa = 10.
dim = n_elements(waka_root)
bslij = fltarr(dim,dimx)
bslky = fltarr(dim,dimx)
vr_md = fltarr(dim,dimx)
vt_md = fltarr(dim,dimx)
vz_md = fltarr(dim,dimx)
pr_md = fltarr(dim,dimx)
rr_md = fltarr(dim,dimx)
mor2 = fltarr(dim,dimx)
mer2 = fltarr(dim,dimx)
R = fltarr(dim,dimx)
sigi = fltarr(dim)
sigk = fltarr(dim)
aa1 = fltarr(dim)
aa0 = fltarr(dim)
aa0[*] = 1.

locout = where(abs(gridx-r0) ge aa,nlocout)
locin = where(abs(gridx-r0) lt aa,nlocin)

nummor2 = (co^2-waka_root^2)*(va^2-waka_root^2)
dnummor2 = (ct^2-waka_root^2)*(co^2+va^2)
nummer2 = (ce^2-waka_root^2)*(vae^2-waka_root^2)
dnummer2 = (cte^2-waka_root^2)*(ce^2+vae^2)
moa2 = ka_root^2*nummor2/dnummor2
mea2 = ka_root^2*nummer2/dnummer2

thrs = 1.e-4
reg0 = where(abs(moa2) lt thrs or abs(mea2) lt thrs)
reg1 = where(moa2 gt thrs and mea2 gt thrs)
reg2 = where(moa2 gt thrs and mea2 lt -thrs)
reg3 = where(moa2 lt -thrs and mea2 gt thrs)
reg4 = where(moa2 lt -thrs and mea2 lt -thrs)

if reg0[0] ne -1 then begin
   aa1[reg0] = 0.
   for i=0,nlocin-1 do begin
      vr_md[reg0,locin[i]] = 0.
      vt_md[reg0,locin[i]] = 1.
      vz_md[reg0,locin[i]] = 0.
      pr_md[reg0,locin[i]] = 0.;ro*aa0/ka_root[reg0]
      rr_md[reg0,locin[i]] = 0.
   endfor
   for i=0,nlocout-1 do begin
      vr_md[reg0,locout[i]] = 0.
      vt_md[reg0,locout[i]] = 1.
      vz_md[reg0,locout[i]] = 0.
      pr_md[reg0,locout[i]] = 0.;re*aa1/ka_root[reg0]
      rr_md[reg0,locout[i]] = 0.
   endfor
endif
if reg1[0] ne -1 then begin
;   dbslij[reg1,*] = dbeseli(sqrt(abs(mor2[reg1,*])),n)
   sigi[reg1] = 1.
   sigk[reg1] = 1.
   aa1[reg1] = aa0[reg1]*beseli(sqrt(abs(moa2[reg1])),n,/double)/beselk(sqrt(abs(mea2[reg1])),n,/double)*(co^2+va^2)*(waka_root[reg1]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg1]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg1,locin[i]] = moa2[reg1]*(abs(gridx[locin[i]]-r0)/aa)^2
;      bslij[reg1,locin[i]] = beseli(sqrt(abs(mor2[reg1,locin[i]])),n,/double)
      vr_md[reg1,locin[i]] = (va^2+co^2)*(waka_root[reg1]^2-ct^2)/(waka_root[reg1]^2*(va^2-waka_root[reg1]^2))*aa0[reg1]*sigi[reg1]*sqrt(abs(moa2[reg1]/aa^2))*dbeseli(sqrt(abs(moa2[reg1]))*abs((gridx[locin[i]]-r0)/aa),n)*(aa/ka_root[reg1])^2
      vt_md[reg1,locin[i]] = (-va^2-co^2+va^2*co^2/waka_root[reg1])*aa0[reg1]*complex(0,1)*n*beseli(sqrt(abs(mor2[reg1,locin[i]])),n)/((waka_root[reg1]*ka_root[reg1])^2-va^2*ka_root[reg1]^2)/abs(gridx[locin[i]]-r0)*aa^2
      vz_md[reg1,locin[i]] = -complex(0,1)/waka_root[reg1]^2/ka_root[reg1]*co^2*aa0[reg1]*beseli(sqrt(abs(mor2[reg1,locin[i]])),n)*aa
      pr_md[reg1,locin[i]] = ro/ka_root[reg1]*(co^2+va^2-va^2*co^2/waka_root[reg1]^2)*beseli(sqrt(abs(mor2[reg1,locin[i]])),n)*aa0[reg1]
      rr_md[reg1,locin[i]] = 1./waka_root[reg1]/ka_root[reg1]*ro*aa0[reg1]*sigi[reg1]*beseli(sqrt(abs(moa2[reg1]))*abs((gridx[locin[i]]-r0)/aa),n,/double)*aa
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg1,locout[i]] = mea2[reg1]*(abs(gridx[locout[i]]-r0)/aa)^2
      vr_md[reg1,locout[i]] = (vae^2+ce^2)*(waka_root[reg1]^2-cte^2)/(waka_root[reg1]^2*(vae^2-waka_root[reg1]^2))*aa1[reg1]*sigk[reg1]*sqrt(abs(mea2[reg1]/aa^2))*dbeselk(sqrt(abs(mea2[reg1]))*abs((gridx[locout[i]]-r0)/aa),n)*(aa/ka_root[reg1])^2
      vt_md[reg1,locout[i]] = (-vae^2-ce^2+vae^2*ce^2/waka_root[reg1]^2)*aa1[reg1]*complex(0,1)*n*beselk(sqrt(abs(mer2[reg1,locout[i]])),n)/((waka_root[reg1]*ka_root[reg1])^2-vae^2*ka_root[reg1]^2)/abs(gridx[locout[i]]-r0)*aa^2
      vz_md[reg1,locout[i]] = -complex(0,1)/waka_root[reg1]^2/ka_root[reg1]*ce^2*aa1[reg1]*beselk(sqrt(abs(mer2[reg1,locout[i]])),n)*aa
      pr_md[reg1,locout[i]] = re/ka_root[reg1]*(ce^2+vae^2-vae^2*ce^2/waka_root[reg1]^2)*beselk(sqrt(abs(mer2[reg1,locout[i]])),n)*aa1[reg1]
      rr_md[reg1,locout[i]] = 1./waka_root[reg1]/ka_root[reg1]*re*aa1[reg1]*sigk[reg1]*beselk(sqrt(abs(mea2[reg1]))*abs((gridx[locout[i]]-r0)/aa),n,/double)*aa
   endfor
endif

if reg2[0] ne -1 then begin
;   dbslij[reg2,*] = dbeseli(sqrt(abs(mor2[reg2,*])),n)
   sigi[reg2] = 1.
   sigk[reg2] = -1.
   aa1[reg2] = aa0[reg2]*beseli(sqrt(abs(moa2[reg2])),n,/double)/besely(sqrt(abs(mea2[reg2])),n,/double)*(co^2+va^2)*(waka_root[reg2]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg2]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg2,locin[i]] = moa2[reg2]*(abs(gridx[locin[i]]-r0)/aa)^2
;      bslij[reg2,locin[i]] = beseli(sqrt(abs(mor2[reg2,locin[i]])),n,/double)
      vr_md[reg2,locin[i]] = (va^2+co^2)*(waka_root[reg2]^2-ct^2)/(waka_root[reg2]^2*(va^2-waka_root[reg2]^2))*aa0[reg2]*sigi[reg2]*sqrt(abs(moa2[reg2]/aa^2))*dbeseli(sqrt(abs(moa2[reg2]))*abs((gridx[locin[i]]-r0)/aa),n)*(aa/ka_root[reg2])^2
      vt_md[reg2,locin[i]] = (-va^2-co^2+va^2*co^2/waka_root[reg2])*aa0[reg2]*complex(0,1)*n*beseli(sqrt(abs(mor2[reg2,locin[i]])),n)/((waka_root[reg2]*ka_root[reg2])^2-va^2*ka_root[reg2]^2)/abs(gridx[locin[i]]-r0)*aa^2
      vz_md[reg2,locin[i]] = -complex(0,1)/waka_root[reg2]^2/ka_root[reg2]*co^2*aa0[reg2]*beseli(sqrt(abs(mor2[reg2,locin[i]])),n)*aa
      pr_md[reg2,locin[i]] = ro/ka_root[reg2]*(co^2+va^2-va^2*co^2/waka_root[reg2]^2)*beseli(sqrt(abs(mor2[reg2,locin[i]])),n)*aa0[reg2]
      rr_md[reg2,locin[i]] = 1./waka_root[reg2]/ka_root[reg2]*ro*aa0[reg2]*sigi[reg2]*beseli(sqrt(abs(moa2[reg2]))*abs((gridx[locin[i]]-r0)/aa),n,/double)*aa
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg2,locout[i]] = mea2[reg2]*(abs(gridx[locout[i]]-r0)/aa)^2
;      bslky[reg2,locout[i]] = besely(sqrt(abs(mer2[reg2,locout[i]])),n,/double)
;      dbslky[reg2,locout[i]] = dbesely(sqrt(abs(mer2[reg2,locout[i]])),n)
      vr_md[reg2,locout[i]] = (vae^2+ce^2)*(waka_root[reg2]^2-cte^2)/(waka_root[reg2]^2*(vae^2-waka_root[reg2]^2))*aa1[reg2]*sigk[reg2]*sqrt(abs(mea2[reg2]/aa^2))*dbesely(sqrt(abs(mea2[reg2]))*abs((gridx[locout[i]]-r0)/aa),n)*(aa/ka_root[reg2])^2
      vt_md[reg2,locout[i]] = (-vae^2-ce^2+vae^2*ce^2/waka_root[reg2]^2)*aa1[reg2]*complex(0,1)*n*besely(sqrt(abs(mer2[reg2,locout[i]])),n)/((waka_root[reg2]*ka_root[reg2])^2-vae^2*ka_root[reg2]^2)/abs(gridx[locout[i]]-r0)*aa^2
      vz_md[reg2,locout[i]] = -complex(0,1)/waka_root[reg2]^2/ka_root[reg2]*ce^2*aa1[reg2]*besely(sqrt(abs(mer2[reg2,locout[i]])),n)*aa
      pr_md[reg2,locout[i]] = re/ka_root[reg2]*(ce^2+vae^2-vae^2*ce^2/waka_root[reg2]^2)*besely(sqrt(abs(mer2[reg2,locout[i]])),n)*aa1[reg2]
      rr_md[reg2,locout[i]] = 1./waka_root[reg2]/ka_root[reg2]*re*aa1[reg2]*sigk[reg2]*besely(sqrt(abs(mea2[reg2]))*abs((gridx[locout[i]]-r0)/aa),n,/double)*aa
   endfor
endif

if reg3[0] ne -1 then begin
;   dbslij[reg3,*] = dbeselj(sqrt(abs(mor2[reg3,*])),n)
   sigi[reg3] = -1.
   sigk[reg3] = 1.
   aa1[reg3] = aa0[reg3]*beselj(sqrt(abs(moa2[reg3])),n,/double)/beselk(sqrt(abs(mea2[reg3])),n,/double)*(co^2+va^2)*(waka_root[reg3]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg3]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg3,locin[i]] = moa2[reg3]*(abs(gridx[locin[i]]-r0)/aa)^2
      R[reg3,locin[i]] = aa0[reg3]*beselj(sqrt(abs(moa2[reg3]))*abs((gridx[locin[i]]-r0)/aa),n,/double)
;      bslij[reg3,locin[i]] = beselj(sqrt(abs(mor2[reg3,locin[i]])),n,/double)
      vr_md[reg3,locin[i]] = (va^2+co^2)*(waka_root[reg3]^2-ct^2)/(waka_root[reg3]^2*(va^2-waka_root[reg3]^2))*aa0[reg3]*sigi[reg3]*sqrt(abs(moa2[reg3]/aa^2))*dbeselj(sqrt(abs(moa2[reg3]))*abs((gridx[locin[i]]-r0)/aa),n)*(aa/ka_root[reg3])^2
      vt_md[reg3,locin[i]] = (-va^2-co^2+va^2*co^2/waka_root[reg3])*aa0[reg3]*complex(0,1)*n*beselj(sqrt(abs(mor2[reg3,locin[i]])),n)/((waka_root[reg3]*ka_root[reg3])^2-va^2*ka_root[reg3]^2)/abs(gridx[locin[i]]-r0)*aa^2
      vz_md[reg3,locin[i]] = -complex(0,1)/waka_root[reg3]^2/ka_root[reg3]*co^2*aa0[reg3]*beselj(sqrt(abs(mor2[reg3,locin[i]])),n)*aa
      pr_md[reg3,locin[i]] = ro/waka_root[reg3]*aa/ka_root[reg3]*(co^2+va^2-va^2*co^2/waka_root[reg3]^2)*beselj(sqrt(abs(moa2[reg3]))*abs((gridx[locin[i]]-r0)/aa),n,/double)*aa0[reg3]
      rr_md[reg3,locin[i]] = 1./waka_root[reg3]*aa/ka_root[reg3]*ro*aa0[reg3]*beselj(sqrt(abs(moa2[reg3]))*abs((gridx[locin[i]]-r0)/aa),n,/double)
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg3,locout[i]] = mea2[reg3]*(abs(gridx[locout[i]]-r0)/aa)^2
      R[reg3,locout[i]] = aa1[reg3]*beselk(sqrt(abs(mea2[reg3]))*abs((gridx[locout[i]]-r0)/aa),n,/double)
;      bslky[reg3,locout[i]] = beselk(sqrt(abs(mer2[reg3,locout[i]])),n,/double)
;      dbslky[reg3,locout[i]] = dbeselk(sqrt(abs(mer2[reg3,locout[i]])),n)
      vr_md[reg3,locout[i]] = (vae^2+ce^2)*(waka_root[reg3]^2-cte^2)/(waka_root[reg3]^2*(vae^2-waka_root[reg3]^2))*aa1[reg3]*sigk[reg3]*sqrt(abs(mea2[reg3]/aa^2))*dbeselk(sqrt(abs(mea2[reg3]))*abs((gridx[locout[i]]-r0)/aa),n)*(aa/ka_root[reg3])^2
      vt_md[reg3,locout[i]] = (-vae^2-ce^2+vae^2*ce^2/waka_root[reg3]^2)*aa1[reg3]*complex(0,1)*n*beselk(sqrt(abs(mer2[reg3,locout[i]])),n)/((waka_root[reg3]*ka_root[reg3])^2-vae^2*ka_root[reg3]^2)/abs(gridx[locout[i]]-r0)*aa^2
      vz_md[reg3,locout[i]] = -complex(0,1)/waka_root[reg3]^2/ka_root[reg3]*ce^2*aa1[reg3]*beselk(sqrt(abs(mer2[reg3,locout[i]])),n)*aa
      pr_md[reg3,locout[i]] = re/waka_root[reg3]*aa/ka_root[reg3]*(ce^2+vae^2-vae^2*ce^2/waka_root[reg3]^2)*beselk(sqrt(abs(mea2[reg3]))*abs((gridx[locout[i]]-r0)/aa),n,/double)*aa1[reg3]
      pr_md[reg3,locout[i]] = re/ka_root[reg3]*(ce^2+vae^2-vae^2*ce^2/waka_root[reg3]^2)*beselk(sqrt(abs(mer2[reg3,locout[i]])),n)*aa1[reg3]
      rr_md[reg3,locout[i]] = 1./waka_root[reg3]*aa/ka_root[reg3]*re*aa1[reg3]*beselk(sqrt(abs(mea2[reg3]))*abs((gridx[locout[i]]-r0)/aa),n,/double)
   endfor
endif

if reg4[0] ne -1 then begin
;   dbslij[reg4,*] = dbeselj(sqrt(abs(mor2[reg4,*])),n)
   sigi[reg4] = -1.
   sigk[reg4] = -1.
   aa1[reg4] = aa0[reg4]*beselj(sqrt(abs(moa2[reg4])),n,/double)/besely(sqrt(abs(mea2[reg4])),n,/double)*(co^2+va^2)*(waka_root[reg4]^2-ct^2)/((ce^2+vae^2)*(waka_root[reg4]^2-cte^2))*ro/re
   for i=0,nlocin-1 do begin
      mor2[reg4,locin[i]] = moa2[reg4]*(abs(gridx[locin[i]]-r0)/aa)^2
 ;     bslij[reg4,locin[i]] = beselj(sqrt(abs(mor2[reg4,locin[i]])),n,/double)
      vr_md[reg4,locin[i]] = (va^2+co^2)*(waka_root[reg4]^2-ct^2)/(waka_root[reg4]^2*(va^2-waka_root[reg4]^2))*aa0[reg4]*sigi[reg4]*sqrt(abs(moa2[reg4]/aa^2))*dbeselj(sqrt(abs(moa2[reg4]))*abs((gridx[locin[i]]-r0)/aa),n)*(aa/ka_root[reg4])^2
      vt_md[reg4,locin[i]] = (-va^2-co^2+va^2*co^2/waka_root[reg4])*aa0[reg4]*complex(0,1)*n*beselj(sqrt(abs(mor2[reg4,locin[i]])),n)/((waka_root[reg4]*ka_root[reg4])^2-va^2*ka_root[reg4]^2)/abs(gridx[locin[i]]-r0)*aa^2
      vz_md[reg4,locin[i]] = -complex(0,1)/waka_root[reg4]^2/ka_root[reg4]*co^2*aa0[reg4]*beselj(sqrt(abs(mor2[reg4,locin[i]])),n)*aa
      pr_md[reg4,locin[i]] = ro/ka_root[reg4]*(co^2+va^2-va^2*co^2/waka_root[reg4]^2)*beselj(sqrt(abs(mor2[reg4,locin[i]])),n)*aa0[reg4]
      rr_md[reg4,locin[i]] = 1./waka_root[reg4]/ka_root[reg4]*ro*aa0[reg4]*sigi[reg4]*beselj(sqrt(abs(moa2[reg4]))*abs((gridx[locin[i]]-r0)/aa),n,/double)*aa
   endfor
   for i=0,nlocout-1 do begin
      mer2[reg4,locout[i]] = mea2[reg4]*(abs(gridx[locout[i]]-r0)/aa)^2
;      bslky[reg4,locout[i]] = besely(sqrt(abs(mer2[reg4,locout[i]])),n,/double)
;      dbslky[reg4,locout[i]] = dbesely(sqrt(abs(mer2[reg4,locout[i]])),n)
      vr_md[reg4,locout[i]] = (vae^2+ce^2)*(waka_root[reg4]^2-cte^2)/(waka_root[reg4]^2*(vae^2-waka_root[reg4]^2))*aa1[reg4]*sigk[reg4]*sqrt(abs(mea2[reg4]/aa^2))*dbesely(sqrt(abs(mea2[reg4]))*abs((gridx[locout[i]]-r0)/aa),n)*(aa/ka_root[reg4])^2
      vt_md[reg4,locout[i]] = (-vae^2-ce^2+vae^2*ce^2/waka_root[reg4]^2)*aa1[reg4]*complex(0,1)*n*besely(sqrt(abs(mer2[reg4,locout[i]])),n)/((waka_root[reg4]*ka_root[reg4])^2-vae^2*ka_root[reg4]^2)/abs(gridx[locout[i]]-r0)*aa^2
      vz_md[reg4,locout[i]] = -complex(0,1)/waka_root[reg4]^2/ka_root[reg4]*ce^2*aa1[reg4]*besely(sqrt(abs(mer2[reg4,locout[i]])),n)*aa
      pr_md[reg4,locout[i]] = re/ka_root[reg4]*(ce^2+vae^2-vae^2*ce^2/waka_root[reg4]^2)*besely(sqrt(abs(mer2[reg4,locout[i]])),n)*aa1[reg4]
      rr_md[reg4,locout[i]] = 1./waka_root[reg4]/ka_root[reg4]*re*aa1[reg4]*sigk[reg4]*besely(sqrt(abs(mea2[reg4]))*abs((gridx[locout[i]]-r0)/aa),n,/double)*aa
   endfor
endif


end
