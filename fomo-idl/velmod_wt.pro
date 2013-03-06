pro velmod_wt, waka_root=waka_root, ka_root=ka_root, gridx=gridx, gridr=gridr, dimt=dimt, reg3=reg3, vr_md=vr_md, aa, wk_rt, ka_rt, kafix, diml, dimz, theta, tarr, gridz, vr_t, vt_t, vz_t, rr_t, rtot_t, ptot_t, dispr, te_t, br_t, bt_t, bz_t, mag=mag, nmode=nmode, no_standing=no_standing,uniform=uniform

if n_params(0) lt 1 then begin
   print,'velmod_wt, waka_root=waka_root, ka_root=ka_root, gridx=gridx, gridr=gridr, dimt=dimt, reg3=reg3, vr_md=vr_md, aa, wk_rt, ka_rt, kafix, diml, dimz, theta, tarr, gridz, vr_t, vt_t, vz_t, rr_t, rtot_t, ptot_t, dispr, te_t, br_t, bt_t, bz_t, mag=mag, nmode=nmode, no_standing=no_standing,uniform=uniform'
   return
endif

; Calculates the modulation in time and space (r,z) of the MHD mode 
; on thermodynamic and geometrical quantities (vr, vt, vz, pr, rr) for
; a specific wavenumber, according to Edwin & Roberts 1983
; INPUT:
; waka_root = w/k roots of dispersion relation
; ka_root = values of ka corresponding to w/k roots (a = radius of cylinder)
; gridx = grid along x axis (perpendicular to cylinder axis)
; dimt = number of points in time dimension
; reg3 = positions in ka_root array corresponding to trapped modes
; vr_md = vr(k,x) : array of radial velocity (wavelength, x position)
; OUTPUT:
; aa = radius of cylinder
; wk_rt = w/k solutions corresponding to trapped modes
; ka_rt = ka values of w/k solutions corresponding to trapped modes
; kafix = location of wavenumber in ka_rt where quantities are calculated
; dimz = dimension in z direction (set equal to 1 wavelength)
; gridz = grid along z axis (longitudinal)
; tarr = time array
; vr_t = vr(dimx,dimz,dimt) : radial velocity
; vz_t = vz(dimx,dimz,dimt) : longitudinal velocity
; rr_t = rr(dimx,dimz,dimt) : density perturbation
; rtot_t = rtot(dimx,dimz,dimt) : total density
; ptot_t = ptot(dimx,dimz,dimt) : total gas pressure
; te_t = te_t(dimx,dimy,dimz) : temperature
; dispr = dispr(dimx,dimz,dimt) : displacement field
; OPTIONAL:
; if mag keyword is set then produces magnetic field cubes:
; br_t = br_t(dimx,dimy,dimz) : radial magnetic field
; bz_t = bz_t(dimx,dimy,dimz) : longitudinal magnetic field
; if keyword 'no_standing' is set then propagating mode case is
; treated, else standing mode is treated

if keyword_set(no_standing) then standing = 0 else standing = 1

common vars1, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, n, ka

if keyword_set(nmode) then n = nmode else n = 0 ; default is sausage mode

co = A1 & va = A2 & ct = A3 & ro = A4 & bo = A5
ce = B1 & vae = B2 & cte = B3 & re = B4 & be = B5

r0 = gridx[n_elements(gridx)-1]/2.
dimx = n_elements(gridx)
dimr = n_elements(gridr)
aa = 10. ; radius of cylinder
if n eq 0 then amp = 0.1 else amp = 0.002 ; amplitude of perturbation
wk_rt = waka_root[reg3] ; reg3 corresponds to the range of coronal solutions
ka_rt = ka_root[reg3]
vr_0 = vr_md[reg3,*]
polyind = 5./3.
po = ro*co^2/polyind
pe = re*ce^2/polyind

;kfixar = where(abs(ka_rt-0.03) lt 0.001) ; for ka closest to 2.244
kfixar = ([min(abs(ka_rt-2.244)),!c])[1]
if kfixar[0] eq -1 then begin print,'no solution for ka' & return & endif
kfloc = ([min(wk_rt[kfixar]),!c])[1]
kafix = kfixar[kfloc]
;kafix = ([min(abs(ka_rt-1.255)),!c])[1] ; for ka closest to 1.255

t_u = 2*!pi/(wk_rt[kafix]*ka_rt[kafix])*aa
tarr = findgen(dimt)/(dimt-1.)*t_u
if n eq 0 then z_u = 2*!pi/(ka_rt[kafix])*aa else z_u = !pi/(ka_rt[kafix])*aa
if keyword_set(uniform) then dimz = round(z_u*dimr/gridr[n_elements(gridr)-1]) else dimz = dimx
gridz = findgen(dimz)/(dimz-1.)*z_u
if n eq 0 then diml = 1 else diml = dimx
l_u = 2*!pi
if n eq 0 then theta = findgen(diml)*l_u else theta = findgen(diml)/diml*l_u
proton = 1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
print,'wavenumber: ',ka_rt[kafix]
print,'phase speed: ',wk_rt[kafix]*100,' km/s'
print,'length of loop: ',z_u/aa,' Mm'
print,'Period of oscillation: ',t_u/60,' min'

stop

; quantities in time:
vr_t = fltarr(dimr,diml,dimz,dimt)
dispr = fltarr(dimr,diml,dimz,dimt)
vt_t = fltarr(dimr,diml,dimz,dimt)
vz_t = fltarr(dimr,diml,dimz,dimt)
rr_t = fltarr(dimr,diml,dimz,dimt)
rtot_t = fltarr(dimr,diml,dimz,dimt)
ptot_t = fltarr(dimr,diml,dimz,dimt)
te_t = fltarr(dimr,diml,dimz,dimt)
if keyword_set(mag) then begin
   br_t = fltarr(dimr,diml,dimz,dimt)
   bt_t = fltarr(dimr,diml,dimz,dimt)
   bz_t = fltarr(dimr,diml,dimz,dimt)
endif

nummor2 = (co^2-wk_rt^2)*(va^2-wk_rt^2)
dnummor2 = (ct^2-wk_rt^2)*(co^2+va^2)
nummer2 = (ce^2-wk_rt^2)*(vae^2-wk_rt^2)
dnummer2 = (cte^2-wk_rt^2)*(ce^2+vae^2)
moa2 = ka_rt^2*nummor2/dnummor2
mea2 = ka_rt^2*nummer2/dnummer2

sigi = -1.
sigk = 1.
;locout = where(abs(gridr) ge aa,nlocout)
;locin = where(abs(gridr) lt aa,nlocin)
;neglocin = where(gridr[locin] lt 0.)
;poslocin = where(gridr[locin] ge 0.)
;neglocout = where(gridr[locout] lt 0.)
;poslocout = where(gridr[locout] ge 0.)
;signr = intarr(dimr)
;signr[locin[poslocin]] = 1
;signr[locin[neglocin]] = -1
;signr[locout[poslocout]] = 1
;signr[locout[neglocout]] = -1
aa0 = 1.

if standing eq 1 then begin
   for k=0,dimt-1 do begin
      for j=0,dimz-1 do begin
         for l=0,diml-1 do begin
            for i=0,dimr-1 do begin 
	       dispr[i,l,j,k] = amp*vr_0[kafix,i]/wk_rt[kafix]*aa/ka_rt[kafix]*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
               displac1 = gridr[i]-dispr[i,l,j,k]
               if (abs(displac1) lt aa) then begin
;               if ((gridr[i] gt aa and abs(displac1) lt aa) or (gridr[i] lt aa and abs(displac1) lt 2*aa)) then begin
                  mor2 = moa2[kafix]*(displac1/aa)^2
                  vr_t[i,l,j,k] = amp*((va^2+co^2)*(wk_rt[kafix]^2-ct^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-va^2))*aa0*sqrt(sigi*moa2[kafix]/aa^2)*dbeselj(sqrt(sigi*mor2),n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  if displac1 eq 0. then begin
                     if n gt 0 then vt_t[i,l,j,k] = amp*(-(va^2+co^2-va^2*co^2/wk_rt[kafix]^2)*aa0*n*1./((wk_rt[kafix]*ka_rt[kafix])^2-va^2*ka_rt[kafix]^2)*aa^2*1./gamma(n+1)*(sqrt(sigi*moa2[kafix]/aa^2)*0.5)^n*(displac1)^(n-1))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*sin(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa) else vt_t[i,l,j,k] = 0.
                  endif else begin
                     vt_t[i,l,j,k] = amp*(-(va^2+co^2-va^2*co^2/wk_rt[kafix]^2)*aa0*n*beselj(sqrt(sigi*mor2),n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-va^2*ka_rt[kafix]^2)/abs(displac1)*aa^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*sin(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  endelse
                  vz_t[i,l,j,k] = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*co^2*aa0*beselj(sqrt(sigi*mor2),n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa)
                  rr_t[i,l,j,k] = amp*(ro/wk_rt[kafix]/ka_rt[kafix]*aa*aa0*beselj(sqrt(sigi*mor2),n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  rtot_t[i,l,j,k] = rr_t[i,l,j,k]+ro
                  ptot_t[i,l,j,k] = amp*(co^2*ro/wk_rt[kafix]/ka_rt[kafix]*aa*aa0*beselj(sqrt(sigi*mor2),n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)+po
                  te_t[i,l,j,k] = ptot_t[i,l,j,k]/rtot_t[i,l,j,k]
                  if keyword_set(mag) then begin
                     br_t[i,l,j,k] = amp*bo/wk_rt[kafix]*((va^2+co^2)*(wk_rt[kafix]^2-ct^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-va^2))*aa0*sqrt(sigi*moa2[kafix]/aa^2)*dbeselj(sqrt(sigi*mor2),n)*(aa/ka_rt[kafix])^2)*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa)
                     if displac1 eq 0. then begin
                        if n gt 0 then bt_t[i,l,j,k] = amp*bo/wk_rt[kafix]*(-(va^2+co^2-va^2*co^2/wk_rt[kafix]^2)*aa0*n*1./((wk_rt[kafix]*ka_rt[kafix])^2-va^2*ka_rt[kafix]^2)*aa^2*1./gamma(n+1)*(sqrt(sigi*moa2[kafix]/aa^2)*0.5)^n*(displac1)^(n-1))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa) else bt_t[i,l,j,k] = 0.
                     endif else begin
                        bt_t[i,l,j,k] = amp*bo/wk_rt[kafix]*(-(va^2+co^2-va^2*co^2/wk_rt[kafix]^2)*aa0*n*beselj(sqrt(sigi*mor2),n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-va^2*ka_rt[kafix]^2)/abs(displac1)*aa^2)*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*sin(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa)
                     endelse
                     bz_t[i,l,j,k] = amp*bo/wk_rt[kafix]/ka_rt[kafix]*aa*(1.-co^2/wk_rt[kafix]^2)*aa0*beselj(sqrt(sigi*mor2),n,/double)*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  endif
               endif else begin
;		  if i ne 0 then if abs(gridr[i-1]-dispr[i-1,l,j,k]) lt aa then loca = i-1
;		  if i ne dimr-1 then if abs(gridr[i+1]-dispr[i+1,l,j,k]) lt aa then loca = i-1
;		  if l eq 0 and j eq 102 and k eq 1 then stop
		  displaca = aa; + dispr[loca,l,j,k]
	          mora2 = moa2[kafix]*(displaca/aa)^2
        	  mera2 = mea2[kafix]*(displaca/aa)^2
;		  aa1 = aa0*dbeselj(sqrt(sigi*mora2),n)/dbeselk(sqrt(sigk*mera2),n)*sigi*moa2[kafix]/(sigk*mea2[kafix])*(wk_rt[kafix]^2-vae^2)/(wk_rt[kafix]^2-va^2)*(co^2+va^2)*(wk_rt[kafix]^2-ct^2)/((ce^2+vae^2)*(wk_rt[kafix]^2-cte^2))
               	  aa1 = aa0*beselj(sqrt(sigi*mora2),n,/double)/beselk(sqrt(sigk*mera2),n,/double)*(co^2+va^2)*(wk_rt[kafix]^2-ct^2)/((ce^2+vae^2)*(wk_rt[kafix]^2-cte^2))*ro/re
                  mer2 = mea2[kafix]*(displac1/aa)^2
                  vr_t[i,l,j,k] = amp*((vae^2+ce^2)*(wk_rt[kafix]^2-cte^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vae^2))*aa1*sqrt(sigk*mea2[kafix]/aa^2)*dbeselk(sqrt(sigk*mer2),n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  vt_t[i,l,j,k] = amp*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[kafix]^2)*aa1*n*beselk(sqrt(sigk*mer2),n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vae^2*ka_rt[kafix]^2)/abs(displac1)*aa^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*sin(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  vz_t[i,l,j,k] = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*ce^2*aa1*beselk(sqrt(sigk*mer2),n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa)
                  rr_t[i,l,j,k] = amp*(re/wk_rt[kafix]/ka_rt[kafix]*aa*aa1*beselk(sqrt(sigk*mer2),n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  rtot_t[i,l,j,k] = rr_t[i,l,j,k]+re
                  ptot_t[i,l,j,k] = amp*(ce^2*re/wk_rt[kafix]/ka_rt[kafix]*aa*aa1*beselk(sqrt(sigk*mer2),n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)+pe
                  te_t[i,l,j,k] = ptot_t[i,l,j,k]/rtot_t[i,l,j,k]
                  if keyword_set(mag) then begin
                     br_t[i,l,j,k] = amp*be/wk_rt[kafix]*((vae^2+ce^2)*(wk_rt[kafix]^2-cte^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vae^2))*aa1*sqrt(sigk*mea2[kafix]/aa^2)*dbeselk(sqrt(sigk*mer2),n)*(aa/ka_rt[kafix])^2)*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa)
                     bt_t[i,l,j,k] = amp*be/wk_rt[kafix]*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[kafix]^2)*aa1*n*beselk(sqrt(sigk*mer2),n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vae^2*ka_rt[kafix]^2)/abs(displac1)*aa^2)*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*sin(n*theta[l])*cos(ka_rt[kafix]*gridz[j]/aa)
                     bz_t[i,l,j,k] = amp*(be/wk_rt[kafix]/ka_rt[kafix]*aa*(1.-ce^2/wk_rt[kafix]^2)*aa1*beselk(sqrt(sigk*mer2),n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(n*theta[l])*sin(ka_rt[kafix]*gridz[j]/aa)
                  endif
               endelse
            endfor
         endfor      
      endfor 
      if n ne 0 then print,string(13b)+' % finished: ',float(k)*100./(dimt-1),format='(a,f4.0,$)'
;stop
   endfor
endif else begin
   for k=0,dimt-1 do begin
      for l=0,diml-1 do begin
         for j=0,dimz-1 do begin
            factora = cos(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa+ka_rt[kafix]*gridz[j]/aa)+complex(0,1)*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa+ka_rt[kafix]*gridz[j]/aa)
            factorb = sin(wk_rt[kafix]*ka_rt[kafix]*tarr[k]/aa)*cos(ka_rt[kafix]*gridz[j]/aa)
            vr_t[i,j,k] = real_part(vr_0[kafix,i]*factora)
            dispr[i,j,k] = vr_t[i,j,k]/wk_rt[kafix]*aa/ka_rt[kafix]*factorb
            if (-dispr[i,j,k]+abs(gridx[i]-r0)) lt aa then begin
               rr_t[i,j,k] = real_part(complex(0,1)/wk_rt[kafix]*aa/ka_rt[kafix]*ro*aa0*beselj(sqrt(abs(moa2[kafix]))*abs((gridx[i]-r0)/aa),n,/double)*factora)
               rtot_t[i,j,k] = real_part(complex(0,1)/wk_rt[kafix]*aa/ka_rt[kafix]*ro*aa0*beselj(sqrt(abs(moa2[kafix]))*abs((gridx[i]-r0)/aa),n,/double)*factora)+ro
               te_t[i,j,k] = rtot_t[i,j,k]^(polyind-1.)
               if keyword_set(mag) then begin
                  br_t[i,j,k] = bo/wk_rt[kafix]*vr_0[kafix,i]*real_part(complex(0,1)*factora)
                  bz_t[i,j,k] = bo/wk_rt[kafix]^3/ka_rt[kafix]*aa*(co^2-wk_rt[kafix]^2)*aa0*beselj(sqrt(abs(moa2[kafix]))*abs((gridx[i]-r0)/aa),n,/double)*real_part(complex(0,1)*factora)
               endif
            endif else begin
               rr_t[i,j,k] = real_part(complex(0,1)/wk_rt[kafix]*aa/ka_rt[kafix]*re*aa1*beselk(sqrt(abs(mea2[kafix]))*abs((gridx[i]-r0)/aa),n,/double)*factora)
               rtot_t[i,j,k] = real_part(complex(0,1)/wk_rt[kafix]*aa/ka_rt[kafix]*re*aa1*beselk(sqrt(abs(mea2[kafix]))*abs((gridx[i]-r0)/aa),n,/double)*factora)+re
               te_t[i,j,k] = rtot_t[i,j,k]^(polyind-1.)
            endelse  
         endfor      
      endfor 
      if n ne 0 then print,string(13b)+' % finished: ',float(k)*100./(dimt-1),format='(a,f4.0,$)'
;stop
   endfor
endelse
vr_t = reform(vr_t)
vt_t = reform(vt_t)
vz_t = reform(vz_t)
rr_t = reform(rr_t)
te_t = reform(te_t)
rtot_t = reform(rtot_t)
ptot_t = reform(ptot_t)
dispr = reform(dispr)
if keyword_set(mag) eq 1 then begin
   br_t = reform(br_t)
   bt_t = reform(bt_t)
   bz_t = reform(bz_t)
endif

end
