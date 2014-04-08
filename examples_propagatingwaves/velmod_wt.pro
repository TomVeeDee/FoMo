pro velmod_wt,waka_root=waka_root,ka_root=ka_root,gridx=gridx,gridr=gridr,dimt=dimt,dimz=dimz,diml=diml,reg3=reg3,vr_md=vr_md,wk_0=wk_0,L=L,aa=aa,wk_rt=wk_rt,ka_rt=ka_rt,kafix=kafix,theta=theta,tarr=tarr,gridz=gridz,mag=mag,nmode=nmode,uniform=uniform,save=save,modelname=modelname,vr_t=vr_t,vt_t=vt_t,vz_t=vz_t,rtot_t=rtot_t,te_t=te_t, btot_t=btot_t

if ~keyword_set(waka_root) then begin
    print,'velmod_wt,waka_root=waka_root,ka_root=ka_root,gridx=gridx,gridr=gridr,dimt=dimt,dimz=dimz,reg3=reg3,vr_md=vr_md,ka_0=ka_0,aa=aa,wk_rt=wk_rt,ka_rt=ka_rt,kafix=kafix,theta=theta,tarr=tarr,gridz=gridz,mag=mag,nmode=nmode,uniform=uniform,save=save,modelname=modelname'
    return
endif

; Calculates the advected modulation in time and space (cylindrical
; coordinates) of a standing MHD mode on thermodynamic and geometrical
; quantities (vr, vt, vz, pr, rr, br, bt, bz) for
; a specific wavenumber, according to Edwin & Roberts 1983
; Given a point in the fixed grid, it finds out by finding the roots
; of equations where the point came from, taking into account
; advection in the radial and azimuthal directions (the longitudinal
; direction is neglected).

; INPUT:
; waka_root = w/k roots of dispersion relation
; ka_root = values of ka corresponding to w/k roots (a = radius of cylinder)
; gridx = grid along x axis (perpendicular to cylinder axis)
; gridr = grid or radial axis
; dimt = number of points in time dimension
; reg3 = positions in ka_root array corresponding to trapped modes
; vr_md = vr(k,x) : array of radial velocity (wavelength, x position)
; wa_0 = desired frequency (omega/kz). Then, the closest wavenumber from wk_root is selected.
; nmode = the azimuthal wavenumber (n=0: sausage mode, n=1: kink mode)

; OUTPUT:
; aa = radius of cylinder
; wk_rt = w/k solutions corresponding to trapped modes
; ka_rt = ka values of w/k solutions corresponding to trapped modes
; kafix = location of wavenumber in ka_rt where quantities are calculated
; dimz = dimension in z direction (set equal to 1 wavelength)
; gridz = grid along z axis (longitudinal)
; theta = grid along the azimuthal direction. It is set to only half
; the cylinder due to axi-symmetry.
; tarr = time array
; If SAVE keyword is set then the following quantities are saved in
; .sav files at each time step:
; vr_t = vr(dimr,diml,dimz) : radial velocity
; vt_t = vr(dimr,diml,dimz) : azimuthal velocity
; vz_t = vz(dimr,diml,dimz) : longitudinal velocity
; rr_t = rr(dimr,diml,dimz) : density perturbation
; rtot_t = rtot(dimr,diml,dimz) : total density
; ptot_t = ptot(dimr,diml,dimz) : total gas pressure
; te_t = te_t(dimr,diml,dimz) : temperature

; OPTIONAL:
; set the keyword UNIFORM if the z axis needs to have the same spatial
; resolution as the x axis (and radial axis). Note that for small
; wavenumbers setting UNIFORM will produce huge arrays in the
; longitudinal direction.
; dimz : dimension of longitudinal axis. Default: dimz = dimx/2
; if mag keyword is set then calculates magnetic field quantities:
; br_t = br_t(dimr,diml,dimz) : radial magnetic field
; bt_t = br_t(dimr,diml,dimz) : azimuthal magnetic field
; bz_t = bz_t(dimr,diml,dimz) : longitudinal magnetic field
; modelname = string with name of model for tagging the output .sav
; files at each time step

if ~keyword_set(modelname) then modnm = '' else modnm = '_'+modelname
  
common vars1, A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, n, ka
common vars2, amp, rad, wk, moaa, meaa, t_k, th_l, z_j, r_i, aa1

if keyword_set(nmode) then n = nmode else n = 0 ; default is sausage mode

co = A1 & vao = A2 & ct = A3 & ro = A4 & bo = A5
ce = B1 & vae = B2 & cte = B3 & re = B4 & be = B5

rad = aa
r0 = gridx[n_elements(gridx)-1]/2.
dimx = n_elements(gridx)
dimr = n_elements(gridr)
if (~(keyword_set(gridy))) then gridy = gridx[dimx/2:*]
dimy = n_elements(gridy)
diml = dimx/2.

if n eq 0 then amp = 2.5 else amp = 5.e-5      ; amplitude of perturbation
wk_rt = waka_root[reg3] ; reg3 corresponds to the range of coronal solutions
ka_rt = ka_root[reg3]
vr_0 = vr_md[reg3,*]
polyind = 5./3.
po = ro*co^2/polyind
pe = re*ce^2/polyind

wfixar = ([min(abs(wk_rt-wk_0)),!c])[1] ; select wa closest to wa_0
if wfixar[0] eq -1 then begin print,'no solution for ka' & return & endif
wfloc = ([min(ka_rt[wfixar]),!c])[1]
wkfix = wfixar[wfloc]

t_u = 2*!pi/(wk_rt[wkfix]*ka_rt[wkfix])*aa      ; Change number of periods here
tarr = findgen(dimt)/(dimt-1.)*t_u
if ~keyword_set(L) then z_u = 15*2*!pi/ka_rt[wkfix] else z_u = L     ; Box dimension in z-direction as a number of wavelengths, if L is not specified

if keyword_set(uniform) and ~keyword_set(dimz) then begin
    dimz = round(z_u*dimx/gridx[n_elements(gridx)-1]) 
endif else begin
    if ~keyword_set(dimz) then dimz = dimx/2
endelse

gridz = findgen(dimz)/(dimz-1.)*z_u

if n eq 0 then diml = 1 else diml = dimx/2.
l_u = !pi
if n eq 0 then theta = findgen(diml)*l_u else theta = findgen(diml)/diml*l_u
wk = wk_rt[wkfix]
ka = ka_rt[wkfix]


proton = 1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
print,'wavenumber: ',ka
print,'phase speed: ',wk*100,' km/s'
print,'vertical wavelength: ',z_u/15,' Mm'
print, 'Loop radius: ', aa,' Mm'
print,'Period of oscillation: ',t_u/2/60,' min'

print,'type .c to continue'
stop

; quantities in time:vr_t = fltarr(dimx,dimy,dimz)
vr_t = fltarr(dimr,diml,dimz,dimt)                    ; For sausage mode we take diml = 1 (later rotate solutions along z-axis)
vt_t = fltarr(dimr,diml,dimz,dimt)                    ; Include diml for consistency with kink mode
vz_t = fltarr(dimr,diml,dimz,dimt)
rr_t = fltarr(dimr,diml,dimz,dimt)
rtot_t = fltarr(dimr,diml,dimz,dimt)
;ptot_t = fltarr(dimr,diml,dimz,dimt)
te_t = fltarr(dimr,diml,dimz,dimt)
if keyword_set(mag) then begin
   br_t = fltarr(dimr,diml,dimz,dimt)
   bt_t = fltarr(dimr,diml,dimz,dimt)
   bz_t = fltarr(dimr,diml,dimz,dimt)
   btot_t = fltarr(dimr,diml,dimz,dimt)
endif

nummor2 = (co^2-wk_rt^2)*(vao^2-wk_rt^2)
dnummor2 = (ct^2-wk_rt^2)*(co^2+vao^2)
nummer2 = (ce^2-wk_rt^2)*(vae^2-wk_rt^2)
dnummer2 = (cte^2-wk_rt^2)*(ce^2+vae^2)
moa2 = ka_rt^2*nummor2/dnummor2
mea2 = ka_rt^2*nummer2/dnummer2

moaa = moa2[wkfix]
meaa = mea2[wkfix]
sigi = -1.
sigk = 1.
aa0 = 1.
counter = 0

for k=0,dimt-1 do begin
    st = string(k,format="(i3.3)")
    for j=0,dimz-1 do begin
        for l=0,diml-1 do begin
            t_k = tarr[k] & th_l = theta[l]  & z_j = gridz[j]
            for i=0,dimr-1 do begin
                r_i = gridr[i]
                X = [r_i,(r_i+aa)/2.,aa]
                if k eq 0 then dr1 = r_i else dr1 = FX_ROOT(X,'disp1',tol=1.e-4,itmax = 10)
                if dr1 ge 0. then begin
                    displac = dr1
                    th_ll = th_l
                endif else begin
                    if k eq 0 then dr2 = r_i else dr2 = FX_ROOT(X,'disp2',tol=1.e-4,itmax = 10)
                    if dr2 eq 0. then begin
                        displac = 0.
                        th_ll = th_l
                    endif else begin
                        if dr2 lt 0. then counter = counter+1
                        displac = dr2
                        th_ll = th_l+!pi 
                    endelse
                endelse
                
                if (abs(displac) lt aa) then begin
                    mor2 = moa2[wkfix]*(displac/aa)^2
                    vr_t[i,l,j,k] = amp*((vao^2+co^2)*(wk_rt[wkfix]^2-ct^2)/(wk_rt[wkfix]^2*(wk_rt[wkfix]^2-vao^2))*aa0*sqrt(sigi*moa2[wkfix]/aa^2)*dbeselj(sqrt(sigi*mor2),n)*(aa/ka_rt[wkfix])^2)*sin(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)          ; Add sin(kz-wt)+sin(kz+wt) to get a standing wave 2*sin(kz)*cos(wt) used in cube with standing wave
                    if displac eq 0. then begin
                        if n gt 0 then vt_t[i,l,j,k] = amp*(-(vao^2+co^2-vao^2*co^2/wk_rt[wkfix]^2)*aa0*n*1./((wk_rt[wkfix]*ka_rt[wkfix])^2-vao^2*ka_rt[wkfix]^2)*aa^2*1./gamma(n+1)*(sqrt(sigi*moa2[wkfix]/aa^2)*0.5)^n*(displac)^(n-1))*cos(wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*sin(n*th_ll)*sin(ka_rt[wkfix]*gridz[j]/aa) else vt_t[i,l,j,k] = 0.
                    endif else begin
                        vt_t[i,l,j,k] = amp*(-(vao^2+co^2-vao^2*co^2/wk_rt[wkfix]^2)*aa0*n*beselj(sqrt(sigi*mor2),n,/double)/((wk_rt[wkfix]*ka_rt[wkfix])^2-vao^2*ka_rt[wkfix]^2)/abs(displac)*aa^2)*cos(wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*sin(n*th_ll)*sin(ka_rt[wkfix]*gridz[j]/aa)
                    endelse
                    vz_t[i,l,j,k] = -amp*(1./wk_rt[wkfix]^2/ka_rt[wkfix]*aa*co^2*aa0*beselj(sqrt(sigi*mor2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)
                    rr_t[i,l,j,k] = amp*(ro/wk_rt[wkfix]/ka_rt[wkfix]*aa*aa0*beselj(sqrt(sigi*mor2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)
                    rtot_t[i,l,j,k] = rr_t[i,l,j,k]+ro
                    ptot_t = amp*(co^2*ro/wk_rt[wkfix]/ka_rt[wkfix]*aa*aa0*beselj(sqrt(sigi*mor2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)+po
                    te_t[i,l,j,k] = ptot_t/rtot_t[i,l,j,k]
                    if keyword_set(mag) then begin
                        br_t[i,l,j,k] = amp*bo/wk_rt[wkfix]*((vao^2+co^2)*(wk_rt[wkfix]^2-ct^2)/(wk_rt[wkfix]^2*(wk_rt[wkfix]^2-vao^2))*aa0*sqrt(sigi*moa2[wkfix]/aa^2)*dbeselj(sqrt(sigi*mor2),n)*(aa/ka_rt[wkfix])^2)*sin(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)
                        if displac eq 0. then begin
                            if n gt 0 then bt_t[i,l,j,k] = amp*bo/wk_rt[wkfix]*(-(vao^2+co^2-vao^2*co^2/wk_rt[wkfix]^2)*aa0*n*1./((wk_rt[wkfix]*ka_rt[wkfix])^2-vao^2*ka_rt[wkfix]^2)*aa^2*1./gamma(n+1)*(sqrt(sigi*moa2[wkfix]/aa^2)*0.5)^n*(displac)^(n-1))*sin(wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)*cos(ka_rt[wkfix]*gridz[j]/aa) else bt_t[i,l,j,k] = 0.
                        endif else begin
                            bt_t[i,l,j,k] = amp*bo/wk_rt[wkfix]*(-(vao^2+co^2-vao^2*co^2/wk_rt[wkfix]^2)*aa0*n*beselj(sqrt(sigi*mor2),n,/double)/((wk_rt[wkfix]*ka_rt[wkfix])^2-vao^2*ka_rt[wkfix]^2)/abs(displac)*aa^2)*sin(wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*sin(n*th_ll)*cos(ka_rt[wkfix]*gridz[j]/aa)
                        endelse
                        bz_t[i,l,j,k] = amp*bo/wk_rt[wkfix]/ka_rt[wkfix]*aa*(1.-co^2/wk_rt[wkfix]^2)*aa0*beselj(sqrt(sigi*mor2),n,/double)*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*th_ll)+bo
                        btot_t[i,l,j,k] = sqrt(br_t[i,l,j,k]^2+bt_t[i,l,j,k]^2+bz_t[i,l,j,k]^2+2*bo*bz_t[i,l,j,k]+bo^2)
                    endif
                endif else begin
                    r_i = displac
                    Y = [th_l,(th_l+!pi/20.),(th_l+!pi/10.)]
                    if displac gt 0. then begin
                        if k eq 0 then dt1 = th_l else dt1 = FX_ROOT(Y,'dispt1',tol=1.e-3,itmax = 10)
                        displact = dt1                
                    endif else begin
                        if displac eq 0. then begin
                            displact = 0.
                        endif else begin
                            if k eq 0 then dt2 = th_l else dt2 = FX_ROOT(Y,'dispt2',tol=1.e-3,itmax = 10)
                            displact = dt2
                        endelse
                    endelse
                    daa = aa
                    mora2 = moa2[wkfix]*(daa/aa)^2
                    mera2 = mea2[wkfix]*(daa/aa)^2
                    aa1 = aa0*beselj(sqrt(sigi*mora2),n,/double)/beselk(sqrt(sigk*mera2),n,/double)*(co^2+vao^2)*(wk_rt[wkfix]^2-ct^2)/((ce^2+vae^2)*(wk_rt[wkfix]^2-cte^2))*ro/re
                    mer2 = mea2[wkfix]*(displac/aa)^2
                    vr_t[i,l,j,k] = amp*((vae^2+ce^2)*(wk_rt[wkfix]^2-cte^2)/(wk_rt[wkfix]^2*(wk_rt[wkfix]^2-vae^2))*aa1*sqrt(sigk*mea2[wkfix]/aa^2)*dbeselk(sqrt(sigk*mer2),n)*(aa/ka_rt[wkfix])^2)*sin(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*displact)
                    vt_t[i,l,j,k] = amp*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[wkfix]^2)*aa1*n*beselk(sqrt(sigk*mer2),n,/double)/((wk_rt[wkfix]*ka_rt[wkfix])^2-vae^2*ka_rt[wkfix]^2)/abs(displac)*aa^2)*sin(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*sin(n*displact)
                    vz_t[i,l,j,k] = -amp*(1./wk_rt[wkfix]^2/ka_rt[wkfix]*aa*ce^2*aa1*beselk(sqrt(sigk*mer2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*displact)
                    rr_t[i,l,j,k] = amp*(re/wk_rt[wkfix]/ka_rt[wkfix]*aa*aa1*beselk(sqrt(sigk*mer2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*displact)
                    rtot_t[i,l,j,k] = rr_t[i,l,j,k]+re
                    ptot_t = amp*(ce^2*re/wk_rt[wkfix]/ka_rt[wkfix]*aa*aa1*beselk(sqrt(sigk*mer2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*displact)+pe
                    te_t[i,l,j,k] = ptot_t/rtot_t[i,l,j,k]
                    if keyword_set(mag) then begin
                        br_t[i,l,j,k] = amp*be/wk_rt[wkfix]*((vae^2+ce^2)*(wk_rt[wkfix]^2-cte^2)/(wk_rt[wkfix]^2*(wk_rt[wkfix]^2-vae^2))*aa1*sqrt(sigk*mea2[wkfix]/aa^2)*dbeselk(sqrt(sigk*mer2),n)*(aa/ka_rt[wkfix])^2)*sin(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*displact)
                        bt_t[i,l,j,k] = amp*be/wk_rt[wkfix]*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[wkfix]^2)*aa1*n*beselk(sqrt(sigk*mer2),n,/double)/((wk_rt[wkfix]*ka_rt[wkfix])^2-vae^2*ka_rt[wkfix]^2)/abs(displac)*aa^2)*sin(wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*sin(n*displact)*cos(ka_rt[wkfix]*gridz[j]/aa)
                        bz_t[i,l,j,k] = amp*(be/wk_rt[wkfix]/ka_rt[wkfix]*aa*(1.-ce^2/wk_rt[wkfix]^2)*aa1*beselk(sqrt(sigk*mer2),n,/double))*cos(ka_rt[wkfix]*gridz[j]/aa-wk_rt[wkfix]*ka_rt[wkfix]*tarr[k]/aa)*cos(n*displact)+be
                        btot_t[i,l,j,k] = sqrt(br_t[i,l,j,k]^2+bt_t[i,l,j,k]^2+bz_t[i,l,j,k]^2+2*be*bz_t[i,l,j,k]+be^2)
                    endif
                endelse
            endfor
        endfor  
      endfor 
      print,string(13b)+' % finished: ',float(k)*100./(dimt-1),format='(a,f4.0,$)'  
 endfor
 if n eq 0 then begin
     vr_t = reform(vr_t[*,*,*,*])   ; Remove the unneccesary l-dimension (we will rotate the solutions along the z-axis and interpolate to obtain
     vt_t = reform(vt_t[*,*,*,*])   ; full data cubes (for sausage mode)
     vz_t = reform(vz_t[*,*,*,*])
     rr_t = reform(rr_t[*,*,*,*])
     rtot_t = reform(rtot_t[*,*,*,*])
;    ptot_t = reform(ptot_t[*,*,*,*])
     te_t = reform(te_t[*,*,*,*])
     if keyword_set(mag) then begin
       br_t = reform(br_t[*,*,*,*])
       bt_t = reform(bt_t[*,*,*,*])
       bz_t = reform(bz_t[*,*,*,*])
       btot_t = reform(btot_t[*,*,*,*])
     endif
endif
if keyword_set(save) then begin
    save,vr_t,vt_t,vz_t,rr_t,rtot_t,te_t,br_t,bt_t,bz_t,btot_t,filename='variables'+modnm+'.sav'
endif else begin
        print,'end of time step. Warning: SAVE keyword not set.'
        stop
endelse
if keyword_set(save) then begin
    save,ro,re,va,vae,co, ce, bo, be, aa,nmode,wk_rt,ka_rt,wkfix,gridx,gridy,gridz,gridr,theta,dimr,dimt,dimz,diml,tarr,nmode,filename='params'+modnm+'.sav'
    print, 'Counter is: '+string(counter)
endif

end
