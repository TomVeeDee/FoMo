pro velmod_wt, waka_root=waka_root, ka_root=ka_root, gridx=gridx, dimt=dimt, reg3=reg3, vr_md=vr_md, aa, wk_rt, ka_rt, k_ind, dimz, tarr, gridz, vr_t, vz_t, rr_t, rtot_t, ptot_t, dispr, te_t, no_standing=no_standing

if n_params(0) lt 1 then begin
   print,'velmod_t, waka_root=waka_root, ka_root=ka_root, gridx=gridx, dimt=dimt, reg3=reg3, vr_md=vr_md, wk_rt, ka_rt, k_ind, dimz, tarr, gridz, vr_t, vz_t, rr_t, rtot_t, ptot_t, dispr, te_t, [no_standing=no_standing]'
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
; k_ind = location of wavenumber in ka_rt where quantities are calculated
; dimz = dimension in z direction (set equal to 1 wavelength)
; gridz = grid along z axis (longitudinal)
; tarr = time array
; vr_t = vr(dimx,dimz,dimt) : radial velocity
; vz_t = vz(dimx,dimz,dimt) : longitudinal velocity
; rr_t = rr(dimx,dimz,dimt) : density perturbation
; rtot_t = rtot(dimx,dimz,dimt) : total density
; ptot_t = ptot(dimx,dimz,dimt) : total pressure
; te_t = te_t(dimx,dimy,dimz) : temperature
; dispr = dispr(dimx,dimz,dimt) : displacement field
; OPTIONAL:
; if keyword 'no_standing' is set then propagating mode case is
; treated, else standing mode is treated

if keyword_set(no_standing) then standing = 0 else standing = 1

common vars1, A1, A2, A3, A4, B1, B2, B3, B4, n, ka

co = A1 & va = A2 & ct = A3 & ro = A4
ce = B1 & vae = B2 & cte = B3 & re = B4

dimx = n_elements(gridx)

n = 0 ; sausage mode
r0 = gridx[-1]/2.
aa = 10. ; radius of cylinder
amp = 0.1 ; amplitude of perturbation
wk_rt = waka_root[reg3]
ka_rt = ka_root[reg3]
vr_0 = vr_md[reg3,*]
gamma = 5./3.
po = ro*co^2/gamma
pe = re*ce^2/gamma

k_ind = ([min(abs(ka_rt-2.244)),!c])[1] ; for ka closest to 2.244
;k_ind = ([min(abs(ka_rt-1.255)),!c])[1] ; for ka closest to 1.255

t_u = 2*!pi/(wk_rt[k_ind]*ka_rt[k_ind])*aa
tarr = findgen(dimt)/(dimt-1.)*t_u
z_u = 2*!pi/(ka_rt[k_ind])*aa
dimz = round(z_u*dimx/gridx[-1])
gridz = findgen(dimz)/(dimz-1.)*z_u
proton = 1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
;mu = 1.27 
;mdivkb = mu*proton/kboltz

; quantities in time:

vr_t = fltarr(dimx,dimz,dimt)
vz_t = fltarr(dimx,dimz,dimt)
rr_t = fltarr(dimx,dimz,dimt)
rtot_t = fltarr(dimx,dimz,dimt)
ptot_t = fltarr(dimx,dimz,dimt)
te_t = fltarr(dimx,dimz,dimt)
dispr = fltarr(dimx,dimz,dimt)
factora = make_array(1,/complex)
factorb = make_array(1,/complex)
aa0 = 1.

nummor2 = (co^2-wk_rt^2)*(va^2-wk_rt^2)
dnummor2 = (ct^2-wk_rt^2)*(co^2+va^2)
nummer2 = (ce^2-wk_rt^2)*(vae^2-wk_rt^2)
dnummer2 = (cte^2-wk_rt^2)*(ce^2+vae^2)
moa2 = ka_rt^2*nummor2/dnummor2
mea2 = ka_rt^2*nummer2/dnummer2

aa1 = aa0*beselj(sqrt(abs(moa2[k_ind])),n,/double)/beselk(sqrt(abs(mea2[k_ind])),n,/double)*(co^2+va^2)*(wk_rt[k_ind]^2-ct^2)/((ce^2+vae^2)*(wk_rt[k_ind]^2-cte^2))*ro/re
for i=0,dimx-1 do begin
   for j=0,dimz-1 do begin
      for k=0,dimt-1 do begin
         if standing eq 1 then begin
            factora = amp*cos(wk_rt[k_ind]*ka_rt[k_ind]*tarr[k]/aa)*(cos(ka_rt[k_ind]*gridz[j]/aa+!pi/2)+complex(0,1)*sin(ka_rt[k_ind]*gridz[j]/aa+!pi/2))
            factorb = amp*sin(wk_rt[k_ind]*ka_rt[k_ind]*tarr[k]/aa)*(cos(ka_rt[k_ind]*gridz[j]/aa+!pi/2)+complex(0,1)*sin(ka_rt[k_ind]*gridz[j]/aa+!pi/2))
            vr_t[i,j,k] = vr_0[k_ind,i]*real_part(factora)
            dispr[i,j,k] = vr_0[k_ind,i]/wk_rt[k_ind]*aa/ka_rt[k_ind]*real_part(factorb)
            if (-dispr[i,j,k]+abs(gridx[i]-r0)) lt aa then begin
               vz_t[i,j,k] = real_part(-complex(0,1)/wk_rt[k_ind]^2*aa/ka_rt[k_ind]*co^2*aa0*beselj(sqrt(abs(moa2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factora)
               rr_t[i,j,k] = -real_part(1./wk_rt[k_ind]*aa/ka_rt[k_ind]*ro*aa0*beselj(sqrt(abs(moa2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factorb)
               rtot_t[i,j,k] = -real_part(1./wk_rt[k_ind]*aa/ka_rt[k_ind]*ro*aa0*beselj(sqrt(abs(moa2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factorb)+ro
               ptot_t[i,j,k] = -real_part(gamma/wk_rt[k_ind]*aa/ka_rt[k_ind]*po*aa0*beselj(sqrt(abs(moa2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factorb)+po
               te_t[i,j,k] = ptot_t[i,j,k]/rtot_t[i,j,k]
;               te_t[i,j,k] = rtot_t[i,j,k]^(gamma-1.)
            endif else begin
               vz_t[i,j,k] = real_part(-complex(0,1)/wk_rt[k_ind]^2*aa/ka_rt[k_ind]*co^2*aa1*beselk(sqrt(abs(mea2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factora)
               rr_t[i,j,k] = -1.*real_part(1./wk_rt[k_ind]*aa/ka_rt[k_ind]*re*aa1*beselk(sqrt(abs(mea2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factorb)
               rtot_t[i,j,k] = -1.*real_part(1./wk_rt[k_ind]*aa/ka_rt[k_ind]*re*aa1*beselk(sqrt(abs(mea2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factorb)+re
               ptot_t[i,j,k] = -1.*real_part(gamma/wk_rt[k_ind]*aa/ka_rt[k_ind]*pe*aa1*beselk(sqrt(abs(mea2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factorb)+pe
               te_t[i,j,k] = ptot_t[i,j,k]/rtot_t[i,j,k]
;               te_t[i,j,k] = rtot_t[i,j,k]^(gamma-1.)
            endelse 
         endif else begin
            factora = cos(wk_rt[k_ind]*ka_rt[k_ind]*tarr[k]/aa+ka_rt[k_ind]*gridz[j]/aa)+complex(0,1)*sin(wk_rt[k_ind]*ka_rt[k_ind]*tarr[k]/aa+ka_rt[k_ind]*gridz[j]/aa)
            factorb = sin(wk_rt[k_ind]*ka_rt[k_ind]*tarr[k]/aa)*cos(ka_rt[k_ind]*gridz[j]/aa)
            vr_t[i,j,k] = real_part(vr_0[k_ind,i]*factora)
            dispr[i,j,k] = vr_t[i,j,k]/wk_rt[k_ind]*aa/ka_rt[k_ind]*factorb
            if (-dispr[i,j,k]+abs(gridx[i]-r0)) lt aa then begin
               rr_t[i,j,k] = real_part(complex(0,1)/wk_rt[k_ind]*aa/ka_rt[k_ind]*ro*aa0*beselj(sqrt(abs(moa2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factora)
               rtot_t[i,j,k] = real_part(complex(0,1)/wk_rt[k_ind]*aa/ka_rt[k_ind]*ro*aa0*beselj(sqrt(abs(moa2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factora)+ro
               te_t[i,j,k] = rtot_t[i,j,k]^(gamma-1.)
            endif else begin
               rr_t[i,j,k] = real_part(complex(0,1)/wk_rt[k_ind]*aa/ka_rt[k_ind]*re*aa1*beselk(sqrt(abs(mea2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factora)
               rtot_t[i,j,k] = real_part(complex(0,1)/wk_rt[k_ind]*aa/ka_rt[k_ind]*re*aa1*beselk(sqrt(abs(mea2[k_ind]))*abs((gridx[i]-r0)/aa),n,/double)*factora)+re
               te_t[i,j,k] = rtot_t[i,j,k]^(gamma-1.)
            endelse  
         endelse             
      endfor
   endfor 
   print,string(13b)+' % finished: ',float(i)*100./(dimx-1),format='(a,f4.0,$)'
endfor

end
