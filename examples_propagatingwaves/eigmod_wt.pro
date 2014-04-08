pro eigmod_wt,waka_root=waka_root,ka_root=ka_root,gridx=gridx,gridr=gridr,dimt=dimt,dimz=dimz,diml=diml,reg3=reg3,vr_md=vr_md, $ 
ka_0=ka_0,aa=aa,wk_rt=wk_rt,ka_rt=ka_rt,kafix=kafix,theta=theta,tarr=tarr,gridz=gridz,mag=mag,nmode=nmode,uniform=uniform,save=save, $
modelname=modelname, vr_t=vr_t,vt_t=vt_t,vz_t=vz_t,rtot_t=rtot_t,te_t=te_t

; Calculates the advected coordinates (Lagrangian frame of reference) of
; the points (plasma elements) inside the loop using all three components 
; velocity. The plasma element coordinates are saved into .dat-files at 
; each time step in order to compute the 3D convex hull in the main function 
; fomo-c/tools of the subdirectory.


; INPUT:
; waka_root = w/k roots of dispersion relation
; ka_root = values of ka corresponding to w/k roots (a = radius of cylinder)
; gridx = grid along x axis (perpendicular to cylinder axis)
; gridr = grid or radial axis
; dimt = number of points in time dimension
; reg3 = positions in ka_root array corresponding to trapped modes
; vr_md = vr(k,x) : array of radial velocity (wavelength, x position)
; ka_0 = desired longitudinal wavenumber. Then, the closest
; wavenumber from ka_root is selected.
; nmode = the azimuthal wavenumber (n=0: sausage mode, n=1: kink mode)

; OUTPUT:
; dimr, dimz, rmax, zmax: fundamental parameters of the fixed grid
; on which we want to calculate the eigenfunctions.
; xpos, ypos, zpos in data files.

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
diml = dimx ; Factor 2 larger since we take the complete cylinder now.

if n eq 0 then amp = 0.1 else amp = 5.e-4 ; amplitude of perturbation
wk_rt = waka_root[reg3] ; reg3 corresponds to the range of coronal solutions
ka_rt = ka_root[reg3]
vr_0 = vr_md[reg3,*]
polyind = 5./3.
po = ro*co^2/polyind
pe = re*ce^2/polyind

kfixar = ([min(abs(ka_rt-ka_0)),!c])[1] ; for ka closest to ka_0
if kfixar[0] eq -1 then begin print,'no solution for ka' & return & endif
kfloc = ([min(wk_rt[kfixar]),!c])[1]
kafix = kfixar[kfloc]

t_u = 2*!pi/(wk_rt[kafix]*ka_rt[kafix])*aa
tarr = findgen(dimt)/(dimt-1.)*t_u
dt = tarr[1] - tarr[0]
if n eq 0 then z_u = !pi/(ka_rt[kafix])*aa else z_u = !pi/(ka_rt[kafix])*aa/2 ; only consider half due to symmetry 

if keyword_set(uniform) and ~keyword_set(dimz) then begin
    dimz = round(z_u*dimx/gridx[n_elements(gridx)-1])
    ;dimz = *dimx/2 
endif else begin
    if ~keyword_set(dimz) then dimz = dimx/2
endelse

gridz = findgen(dimz)/(dimz-1.)*z_u

if n eq 0 then diml = 1 else diml = dimx
l_u = 2*!pi ; only half the cylinder is necessary due to symmetry, but now we try the complete cylinder
if n eq 0 then theta = findgen(diml)*l_u else theta = findgen(diml)/diml*l_u
wk = wk_rt[kafix]
ka = ka_rt[kafix]

nummor2 = (co^2-wk_rt^2)*(vao^2-wk_rt^2)
dnummor2 = (ct^2-wk_rt^2)*(co^2+vao^2)
nummer2 = (ce^2-wk_rt^2)*(vae^2-wk_rt^2)
dnummer2 = (cte^2-wk_rt^2)*(ce^2+vae^2)
moa2 = ka_rt^2*nummor2/dnummor2
mea2 = ka_rt^2*nummer2/dnummer2

moaa = moa2[kafix]  ;mora2 is the same
meaa = mea2[kafix]  ;mera2 is the same
sigi = -1.
sigk = 1.
aa0 = 1.
aa1 = aa0*beselj(sqrt(sigi*moaa),n,/double)/beselk(sqrt(sigk*meaa),n,/double)*(co^2+vao^2)*(wk_rt[kafix]^2-ct^2)/((ce^2+vae^2)*(wk_rt[kafix]^2-cte^2))*ro/re

proton = 1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
print,'wavenumber: ',ka
print,'phase speed: ',wk*100,' km/s'
print,'length of loop: ',z_u/aa,' Mm'
print,'Period of oscillation: ',t_u/60,' min'

; test problem with small values of dimr, diml, dimz
;dimr = 10
;diml = 9
;dimz = 7
dimt = 43

vri_t = fltarr(dimt, dimr*diml*dimz)
vti_t = fltarr(dimt, dimr*diml*dimz)
vzi_t = fltarr(dimt, dimr*diml*dimz)
vre_t = fltarr(dimt, dimr*diml*dimz)
vte_t = fltarr(dimt, dimr*diml*dimz)
vze_t = fltarr(dimt, dimr*diml*dimz)
loc = fltarr(dimr*diml*dimz)

rpos = fltarr(dimt, dimr*diml*dimz)
thpos = fltarr(dimt, dimr*diml*dimz)
zpos = fltarr(dimt, dimr*diml*dimz)

; Write information concerning the grid to a small file, used in c++ to determine
; the Eulerian grid points to write the eigenfunctions into.

OPENW, lun, '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/alldata.dat', /GET_LUN
printf, lun, co
printf, lun, vao
printf, lun, ct
printf, lun, ro
printf, lun, bo
printf, lun, po
printf, lun, ce
printf, lun, vae
printf, lun, cte
printf, lun, re
printf, lun, be
printf, lun, pe
printf, lun, aa
printf, lun, amp
printf, lun, kafix
printf, lun, ka_rt[kafix]
printf, lun, wk_rt[kafix]
printf, lun, moaa
printf, lun, meaa
printf, lun, aa1
printf, lun, t_u
printf, lun, dt
printf, lun, dimr
printf, lun, diml
printf, lun, dimz
printf, lun, dimt
CLOSE, lun
FREE_LUN, lun

print, 'Stop here for C++ colour interpolation'
print, 'To advect using IDL, continue'
stop

; Intermediate step: we quit the CGAL paths and just interpolate the position as a colour.
; Below is interpolation using IDL.



















OPENW, lun, '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/griddata.dat', /GET_LUN
printf, lun, dimr
printf, lun, dimz
printf, lun, gridr[dimr-1]
printf, lun, gridz[dimz-1]
CLOSE, lun
FREE_LUN, lun

; Convert regular grid into coordinates of irregularly spaced points (for use in fomo-c after advection)
; Compute initial values of perturbed thermodynamic quantities and velocity
; Save only the point coordinates, converted into Cartesian, in a file, since the other TD quantities have been checked in test runs
; We only need the interior points to determine which grid points are inside or outside the loop

 filename = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoord'+string(0,format="(i3.3)")+'.dat'
 filenamel = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoordleft'+string(0,format="(i3.3)")+'.dat'
 filenamer = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoordright'+string(0,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
 OPENW, lunl, filenamel, /GET_LUN
 OPENW, lunr, filenamer, /GET_LUN

; Initialisation of the point coordinates (the points to be advected) from the grid and subdivision into left, loop, and right

counter = 0.
counterl = 0.
counterr = 0.

for i=0, dimr-1 do begin
  for j = 0, dimz-1 do begin
    for k=0, diml-1 do begin
        rpos[0, dimz*diml*i + diml*j + k] = gridr[i]
        thpos[0, dimz*diml*i + diml*j + k] = theta[k]
        zpos[0, dimz*diml*i + diml*j + k] = gridz[j]     
      if rpos[0, dimz*diml*i + diml*j + k] lt aa then begin  
        printf, lun, rpos[0, dimz*diml*i + diml*j + k]*cos(thpos[0, dimz*diml*i + diml*j + k])
        printf, lun, rpos[0, dimz*diml*i + diml*j + k]*sin(thpos[0, dimz*diml*i + diml*j + k])
        printf, lun, zpos[0, dimz*diml*i + diml*j + k]
        loc[dimz*diml*i + diml*j + k] = 0
        counter = counter + 1
      endif else begin 
        if !pi/2 lt thpos[0, dimz*diml*i + diml*j + k] and thpos[0, dimz*diml*i + diml*j + k] lt 3*!pi/2 then begin
        printf, lunl, rpos[0, dimz*diml*i + diml*j + k]*cos(thpos[0, dimz*diml*i + diml*j + k])
        printf, lunl, rpos[0, dimz*diml*i + diml*j + k]*sin(thpos[0, dimz*diml*i + diml*j + k])
        printf, lunl, zpos[0, dimz*diml*i + diml*j + k]
        counterl = counterl + 1
        loc[dimz*diml*i + diml*j + k] = -1
      endif else begin
        printf, lunr, rpos[0, dimz*diml*i + diml*j + k]*cos(thpos[0, dimz*diml*i + diml*j + k])
        printf, lunr, rpos[0, dimz*diml*i + diml*j + k]*sin(thpos[0, dimz*diml*i + diml*j + k])
        printf, lunr, zpos[0, dimz*diml*i + diml*j + k]
        counterr = counterr + 1
        loc[dimz*diml*i + diml*j + k] = 1   
        endelse
      endelse
    endfor
  endfor       
endfor

CLOSE, lun
FREE_LUN, lun
CLOSE, lunl
FREE_LUN, lunl
CLOSE, lunr
FREE_LUN, lunr

; Calculate velocity, advect points, write out points for next time steps.

for t = 1, dimt-1. do begin
  
 filenamet = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoord'+string(t,format="(i3.3)")+'.dat'
 filenametl = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoordleft'+string(t,format="(i3.3)")+'.dat'
 filenametr = '/users/cpa/sgijsen/fomo/version_patrick_nov13/examples/data/advectedeigf/largeba5104completecyl/advectedcoordright'+string(t,format="(i3.3)")+'.dat'
 OPENW, lun, filenamet, /GET_LUN
 OPENW, lunl, filenametl, /GET_LUN
 OPENW, lunr, filenametr, /GET_LUN

for i=0, dimr-1 do begin
  for j = 0, dimz-1 do begin
    for k=0, diml-1 do begin
  
  if loc[dimz*diml*i + diml*j + k] eq 0 then begin
    vri_t[t-1, dimz*diml*i + diml*j + k] = amp*((vao^2+co^2)*(wk_rt[kafix]^2-ct^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vao^2))*aa0*sqrt(sigi*moaa/aa^2)*dbeselj(sqrt(sigi*moaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*cos(n*thpos[t-1, dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    vti_t[t-1, dimz*diml*i + diml*j + k] = amp*(-(vao^2+co^2-vao^2*co^2/wk_rt[kafix]^2)*aa0*n*beselj(sqrt(sigi*moaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vao^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*sin(n*thpos[t-1, dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    vzi_t[t-1, dimz*diml*i + diml*j + k] = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*co^2*aa0*beselj(sqrt(sigi*moaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*cos(n*thpos[t-1, dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    rpos[t, dimz*diml*i + diml*j + k] = rpos[t-1, dimz*diml*i + diml*j + k] + vri_t[t-1, dimz*diml*i + diml*j + k]*dt
    thpos[t, dimz*diml*i + diml*j + k] = thpos[t-1, dimz*diml*i + diml*j + k] + vti_t[t-1, dimz*diml*i + diml*j + k]*dt
    zpos[t,dimz*diml*i + diml*j + k] = zpos[t-1, dimz*diml*i + diml*j + k] + vzi_t[t-1, dimz*diml*i + diml*j + k]*dt  
    printf, lun, rpos[t, dimz*diml*i + diml*j + k]*cos(thpos[t, dimz*diml*i + diml*j + k])
    printf, lun, rpos[t, dimz*diml*i + diml*j + k]*sin(thpos[t, dimz*diml*i + diml*j + k])
    printf, lun, zpos[t, dimz*diml*i + diml*j + k]
  endif
  
  if loc[dimz*diml*i + diml*j + k] eq -1 then begin
    vre_t[t-1, dimz*diml*i + diml*j + k] = amp*((vae^2+ce^2)*(wk_rt[kafix]^2-cte^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vae^2))*aa1*sqrt(sigk*meaa/aa^2)*dbeselk(sqrt(sigk*meaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*cos(n*thpos[t-1, dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    vte_t[t-1, dimz*diml*i + diml*j + k] = amp*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[kafix]^2)*aa1*n*beselk(sqrt(sigk*meaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vae^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*sin(n*thpos[t-1, dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    vze_t[t-1, dimz*diml*i + diml*j + k] = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*ce^2*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*cos(n*thpos[t-1, dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    rpos[t, dimz*diml*i + diml*j + k] = rpos[t-1, dimz*diml*i + diml*j + k] + vri_t[t-1, dimz*diml*i + diml*j + k]*dt
    thpos[t, dimz*diml*i + diml*j + k] = thpos[t-1, dimz*diml*i + diml*j + k] + vti_t[t-1, dimz*diml*i + diml*j + k]*dt
    zpos[t,dimz*diml*i + diml*j + k] = zpos[t-1, dimz*diml*i + diml*j + k] + vzi_t[t-1, dimz*diml*i + diml*j + k]*dt  
    printf, lunl, rpos[t, dimz*diml*i + diml*j + k]*cos(thpos[t, dimz*diml*i + diml*j + k])
    printf, lunl, rpos[t, dimz*diml*i + diml*j + k]*sin(thpos[t, dimz*diml*i + diml*j + k])
    printf, lunl, zpos[t, dimz*diml*i + diml*j + k]    
  endif
    
  if loc[dimz*diml*i + diml*j + k] eq 1 then begin
    vre_t[t-1, dimz*diml*i + diml*j + k] = amp*((vae^2+ce^2)*(wk_rt[kafix]^2-cte^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vae^2))*aa1*sqrt(sigk*meaa/aa^2)*dbeselk(sqrt(sigk*meaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*cos(n*thpos[t-1, dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    vte_t[t-1, dimz*diml*i + diml*j + k] = amp*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[kafix]^2)*aa1*n*beselk(sqrt(sigk*meaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vae^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*sin(n*thpos[t-1, dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    vze_t[t-1, dimz*diml*i + diml*j + k] = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*ce^2*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[t-1, dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t-1]/aa)*cos(n*thpos[t-1, dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[t-1, dimz*diml*i + diml*j + k]/aa)
    rpos[t, dimz*diml*i + diml*j + k] = rpos[t-1, dimz*diml*i + diml*j + k] + vri_t[t-1, dimz*diml*i + diml*j + k]*dt
    thpos[t, dimz*diml*i + diml*j + k] = thpos[t-1, dimz*diml*i + diml*j + k] + vti_t[t-1, dimz*diml*i + diml*j + k]*dt
    zpos[t,dimz*diml*i + diml*j + k] = zpos[t-1, dimz*diml*i + diml*j + k] + vzi_t[t-1, dimz*diml*i + diml*j + k]*dt  
    printf, lunr, rpos[t, dimz*diml*i + diml*j + k]*cos(thpos[t, dimz*diml*i + diml*j + k])
    printf, lunr, rpos[t, dimz*diml*i + diml*j + k]*sin(thpos[t, dimz*diml*i + diml*j + k])
    printf, lunr, zpos[t, dimz*diml*i + diml*j + k]    
  endif  
endfor
endfor
endfor      
  print, 'time step'+strtrim(string(t),1) 
  CLOSE, lun
  FREE_LUN, lun
  CLOSE, lunl
  FREE_LUN, lunl
  CLOSE, lunr
  FREE_LUN, lunr
endfor
end