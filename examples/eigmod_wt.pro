pro eigmod_wt,waka_root=waka_root,ka_root=ka_root,gridx=gridx,gridr=gridr,dimx=dimx,dimt=dimt,dimz=dimz,diml=diml,reg3=reg3,vr_md=vr_md, $ 
ka_0=ka_0,aa=aa,wk_rt=wk_rt,ka_rt=ka_rt,kafix=kafix,theta=theta,tarr=tarr,gridz=gridz,mag=mag,nmode=nmode,uniform=uniform,save=save, $
modelname=modelname, vr_t=vr_t,vt_t=vt_t,vz_t=vz_t,rtot_t=rtot_t,te_t=te_t

; Calculates the advected coordinates (Lagrangian frame of reference) of
; the points (plasma elements) inside the loop using all three components 
; velocity. The plasma element coordinates are saved into .dat-files at 
; each time step in order to compute the 3D convex hull in the main function 
; fomo-c/tools of the subdirectory. (Situation Feb 14)

; April 14: We calculate all necessary eigenfunctions in Lagrangian frame of reference
; (Irregular grid)


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

if n eq 0 then amp = 0.1 else amp = 5.e-1 ; amplitude of perturbation
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

if ~keyword_set(L) then z_u = 2*!pi/ka_rt[kafix]*aa else z_u = L*aa     ; Box dimension in z-direction as a number of wavelengths (units of 100km), if L is not specified
gridz = findgen(dimz)/(dimz-1.)*z_u

rad = aa
gridx = findgen(dimx)*gridz[1]
r0 = gridx[n_elements(gridx)-1]/2.
dimr = dimx
gridr = findgen(dimr)/(dimr-1)*abs(r0-gridx[dimx-1])
if (~(keyword_set(gridy))) then gridy = gridx[dimx/2:*]
dimy = n_elements(gridy)
diml = dimx ; Factor 2 larger since we take the complete cylinder now (could be changed in future).

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
norm = 1.e5
cgsfactor = 1.e6
factor = norm^2*proton*cgsfactor

print,'wavenumber: ',ka
print,'phase speed: ',wk*100,' km/s'
print,'Height of box (z): ',z_u/10.,' Mm'
print, 'Loop radius: ', aa/10,' Mm'
print,'Period of oscillation: ',t_u/2/60,' min'
print, 'Resolution: xres = yres = ',gridx[1]/10, ' Mm'
print, 'zres = ', gridz[1]/10, ' Mm'
print, 'To continue, type .c'
stop

; test problem with small values of dimr, diml, dimz
;dimr = 100
;diml = 9
;dimz = 7
;dimt = 7

;vri_t = fltarr(dimr*diml*dimz)
;vti_t = fltarr(dimr*diml*dimz)
;vzi_t = fltarr(dimr*diml*dimz)
;vre_t = fltarr(dimr*diml*dimz)
;vte_t = fltarr(dimr*diml*dimz)
;vze_t = fltarr(dimr*diml*dimz)
loc = fltarr(dimr*diml*dimz)

rpos = fltarr(dimr*diml*dimz)
thpos = fltarr(dimr*diml*dimz)
zpos = fltarr(dimr*diml*dimz)



; Write information concerning the grid to a small file to determine the Eulerian grid points to write the eigenfunctions into.

OPENW, lun, 'advectedeigfA05/initialdata.dat', /GET_LUN
printf, lun, 'Interior sound speed, alfven, cusp, density, mag, pressure'
printf, lun, 'Same for outside'
printf, lun, 'R, amplitude, wavenumber index, phase speed, frequency'
printf, lun, 'kappa inside and outside, aa1, period, time step'
printf, lun, 'dimr diml dimz dimt'
printf, lun, co, vao, ct, ro, bo, po
printf, lun, ce, vae, cte, re, be, pe
printf, lun, aa, amp, kafix, ka_rt[kafix], wk_rt[kafix]
printf, lun, moaa, meaa, aa1, t_u, dt
printf, lun, dimr, diml, dimz, dimt
CLOSE, lun
FREE_LUN, lun





; Convert regular grid into coordinates of irregularly spaced points (for use in fomo-c after advection)
; Compute initial values of perturbed thermodynamic quantities and velocity

; Time zero


 filename = '/users/cpa/sgijsen/fomo/stiefApr1614/examples/advectedeigfA05/eigft'+string(0,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
 printf, lun, 3
 printf, lun, 8812
 printf, lun, dimx*dimy*dimz
 printf, lun, 5

; Initialisation of the point coordinates (the points to be advected) from the grid and subdivision into left, loop, and right

countero = 0.
countere = 0.

for i=0, dimr-1 do begin
  for j = 0, dimz-1 do begin
    for k=0, diml-1 do begin
        rpos[dimz*diml*i + diml*j + k] = gridr[i]
        thpos[dimz*diml*i + diml*j + k] = theta[k]
        zpos[dimz*diml*i + diml*j + k] = gridz[j]     
        xpos = rpos[dimz*diml*i + diml*j + k]*cos(thpos[dimz*diml*i + diml*j + k])
        ypos = rpos[dimz*diml*i + diml*j + k]*sin(thpos[dimz*diml*i + diml*j + k])
        ;zposq = zpos[dimz*diml*i + diml*j + k]
        if rpos[dimz*diml*i + diml*j + k] lt aa then begin  ; specify loc, calculate eigenf, save position and eigenf, advect
          loc[dimz*diml*i + diml*j + k] = 0
          countero = countero + 1
          
          vri_t = amp*((vao^2+co^2)*(wk_rt[kafix]^2-ct^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vao^2))*aa0*sqrt(sigi*moaa/aa^2)*dbeselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vti_t = amp*(-(vao^2+co^2-vao^2*co^2/wk_rt[kafix]^2)*aa0*n*beselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vao^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*sin(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vzi_t = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*co^2*aa0*beselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vxi_t = vri_t * cos(thpos[dimz*diml*i+diml*j+k]) - vti_t*rpos[dimz*diml*i+diml*j+k]*sin(thpos[dimz*diml*i+diml*j+k])
          vyi_t = vri_t * sin(thpos[dimz*diml*i+diml*j+k]) + vti_t*rpos[dimz*diml*i+diml*j+k]*cos(thpos[dimz*diml*i+diml*j+k])
          rhi_t = amp*(ro/wk_rt[kafix]/ka_rt[kafix]*aa*aa0*beselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)+ro
          pri_t = amp*(co^2*ro/wk_rt[kafix]/ka_rt[kafix]*aa*aa0*beselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)+po
          tei_t = pri_t/rhi_t
          
          printf, lun, xpos, ypos, zpos[dimz*diml*i + diml*j + k], rhi_t/factor, tei_t*norm, vxi_t, vyi_t, vzi_t
          
          rpos[dimz*diml*i+diml*j+k] = rpos[dimz*diml*i+diml*j+k] + vri_t*dt
          thpos[dimz*diml*i+diml*j+k] = thpos[dimz*diml*i+diml*j+k] + vti_t*dt
          zpos[dimz*diml*i+diml*j+k] = zpos[dimz*diml*i+diml*j+k] + vzi_t*dt
          
        endif else begin
          loc[dimz*diml*i + diml*j + k] = 1
          countere = countere + 1
          
          vre_t = amp*((vae^2+ce^2)*(wk_rt[kafix]^2-cte^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vae^2))*aa1*sqrt(sigk*meaa/aa^2)*dbeselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vte_t = amp*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[kafix]^2)*aa1*n*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vae^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*sin(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vze_t = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*ce^2*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vxe_t = vre_t * cos(thpos[dimz*diml*i+diml*j+k]) - vte_t*rpos[dimz*diml*i+diml*j+k]*sin(thpos[dimz*diml*i+diml*j+k])
          vye_t = vre_t * sin(thpos[dimz*diml*i+diml*j+k]) + vte_t*rpos[dimz*diml*i+diml*j+k]*cos(thpos[dimz*diml*i+diml*j+k])
          rhe_t = amp*(re/wk_rt[kafix]/ka_rt[kafix]*aa*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)+re
          pre_t = amp*(ce^2*re/wk_rt[kafix]/ka_rt[kafix]*aa*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)+pe
          tee_t = pre_t/rhe_t
          
          printf, lun, xpos, ypos, zpos[dimz*diml*i + diml*j + k], rhe_t/factor, tee_t*norm, vxe_t, vye_t, vze_t
          
          rpos[dimz*diml*i+diml*j+k] = rpos[dimz*diml*i+diml*j+k] + vre_t*dt
          thpos[dimz*diml*i+diml*j+k] = thpos[dimz*diml*i+diml*j+k] + vte_t*dt
          zpos[dimz*diml*i+diml*j+k] = zpos[dimz*diml*i+diml*j+k] + vze_t*dt          
      endelse
    endfor
  endfor       
endfor

CLOSE, lun
FREE_LUN, lun

print, 'Number of points inside' + string(countero)
print, 'Number of points outside' + string(countere)

; Calculate velocity [eigenfunctions], advect points, write out points for next time steps.

for t = 1, dimt-1. do begin
  
 filename = '/users/cpa/sgijsen/fomo/stiefApr1614/examples/advectedeigfA05/eigft'+string(t,format="(i3.3)")+'.dat'
 OPENW, lun, filename, /GET_LUN
 printf, lun, 3
 printf, lun, 8812+t
 printf, lun, dimx*dimy*dimz
 printf, lun, 5

for i=0, dimr-1 do begin
  for j = 0, dimz-1 do begin
    for k=0, diml-1 do begin
      xpos = rpos[dimz*diml*i + diml*j + k]*cos(thpos[dimz*diml*i + diml*j + k])
      ypos = rpos[dimz*diml*i + diml*j + k]*sin(thpos[dimz*diml*i + diml*j + k])  
      if loc[dimz*diml*i + diml*j + k] eq 0 then begin ; Inside the loop
          vri_t = amp*((vao^2+co^2)*(wk_rt[kafix]^2-ct^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vao^2))*aa0*sqrt(sigi*moaa/aa^2)*dbeselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vti_t = amp*(-(vao^2+co^2-vao^2*co^2/wk_rt[kafix]^2)*aa0*n*beselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vao^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t]/aa)*sin(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vzi_t = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*co^2*aa0*beselj(sqrt(sigi*moaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[t]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vxi_t = vri_t * cos(thpos[dimz*diml*i+diml*j+k]) - vti_t*rpos[dimz*diml*i+diml*j+k]*sin(thpos[dimz*diml*i+diml*j+k])
          vyi_t = vri_t * sin(thpos[dimz*diml*i+diml*j+k]) + vti_t*rpos[dimz*diml*i+diml*j+k]*cos(thpos[dimz*diml*i+diml*j+k])
          
          printf, lun, xpos, ypos, zpos[dimz*diml*i + diml*j + k], rhi_t/factor, tei_t*norm, vxi_t, vyi_t, vzi_t
          
          rpos[dimz*diml*i+diml*j+k] = rpos[dimz*diml*i+diml*j+k] + vri_t*dt
          thpos[dimz*diml*i+diml*j+k] = thpos[dimz*diml*i+diml*j+k] + vti_t*dt
          zpos[dimz*diml*i+diml*j+k] = zpos[dimz*diml*i+diml*j+k] + vzi_t*dt
      endif
    
      if loc[dimz*diml*i + diml*j + k] eq 1 then begin ; Outside the loop
          vre_t = amp*((vae^2+ce^2)*(wk_rt[kafix]^2-cte^2)/(wk_rt[kafix]^2*(wk_rt[kafix]^2-vae^2))*aa1*sqrt(sigk*meaa/aa^2)*dbeselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n)*(aa/ka_rt[kafix])^2)*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vte_t = amp*(-(vae^2+ce^2-vae^2*ce^2/wk_rt[kafix]^2)*aa1*n*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double)/((wk_rt[kafix]*ka_rt[kafix])^2-vae^2*ka_rt[kafix]^2))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*sin(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vze_t = amp*(1./wk_rt[kafix]^2/ka_rt[kafix]*aa*ce^2*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*cos(wk_rt[kafix]*ka_rt[kafix]*tarr[0]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*cos(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)
          vxe_t = vre_t * cos(thpos[dimz*diml*i+diml*j+k]) - vte_t*rpos[dimz*diml*i+diml*j+k]*sin(thpos[dimz*diml*i+diml*j+k])
          vye_t = vre_t * sin(thpos[dimz*diml*i+diml*j+k]) + vte_t*rpos[dimz*diml*i+diml*j+k]*cos(thpos[dimz*diml*i+diml*j+k])
          rhe_t = amp*(re/wk_rt[kafix]/ka_rt[kafix]*aa*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[t]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)+re
          pre_t = amp*(ce^2*re/wk_rt[kafix]/ka_rt[kafix]*aa*aa1*beselk(sqrt(sigk*meaa/aa^2)*rpos[dimz*diml*i + diml*j + k],n,/double))*sin(wk_rt[kafix]*ka_rt[kafix]*tarr[t]/aa)*cos(n*thpos[dimz*diml*i + diml*j + k])*sin(ka_rt[kafix]*zpos[dimz*diml*i + diml*j + k]/aa)+pe
          tee_t = pre_t/rhe_t
          
          printf, lun, xpos, ypos, zpos[dimz*diml*i + diml*j + k], rhe_t/factor, tee_t*norm, vxe_t, vye_t, vze_t
          
          rpos[dimz*diml*i+diml*j+k] = rpos[dimz*diml*i+diml*j+k] + vre_t*dt
          thpos[dimz*diml*i+diml*j+k] = thpos[dimz*diml*i+diml*j+k] + vte_t*dt
          zpos[dimz*diml*i+diml*j+k] = zpos[dimz*diml*i+diml*j+k] + vze_t*dt   
      endif  
    endfor
  endfor
  if i MOD 20 eq 0 then print, 'step '+string(i)+' of r-position'
endfor      
  print, 'time step'+strtrim(string(t),1) 
  CLOSE, lun
  FREE_LUN, lun
endfor

print, 'Eigenfunctions written to file for fomo-c.'
end