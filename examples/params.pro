
pro params, set=set, it=it, ix=ix, ro=ro, re=re, valfven=va, vae=vae, co=co, ce=ce, aa=aa, r0=r0, gridx=gridx, gridy=gridy, gridz=gridz, dimx=dimx, dimy=dimy, dimz=dimz, tarr=tarr, kafix=kafix, ka_rt=ka_rt, wk_rt=wk_rt, te=te, rho=rho, vr=vr, vz=vz, velx=velx, vely=vely, velz=velz

if keyword_set(set) eq 0 then begin
   print,'Check input and output directories first'
   print,'params, set=set, it=it, ix=ix, ro=ro, re=re, valfven=va, vae=vae, co=co, ce=ce, aa=aa, r0=r0, gridx=gridx, gridy=gridy, gridz=gridz, dimx=dimx, dimy=dimy, dimz=dimz, tarr=tarr, kafix=kafix, ka_rt=ka_rt, wk_rt=wk_rt, te=te, rho=rho, vr=vr, vz=vz, velx=velx, vely=vely, velz=velz'
   return
endif

; Read paramaters from params_'+kanm+'.sav file and from one
; slice'...'.sav file determined by set, ix and it
; INPUT:
; set = case considered: 
;       set = 2 -> base model ka = 2.24
;       set = 21 -> high T model (base model with high external temp.)
;       set = 22 -> smooth model (base model with smooth density prof.)
;       set = 24 -> hgres model (base model with high spatial resolution)
;       set = 3 -> long lambda model (longer wavelength, ka = 1.25)
; it = time step
; ix = position along x axis (radial axis)
; OUTPUT:
; ro = internal density: rho[dimx/2,dimy/2,dimz/2]
; re = external density: rho[0,0,0]
; valfven = internal Alfven velocity: calculated at (dimx/2,dimy/2,dimz/2)
; vae = external Alfven velocity: calculated at (0,0,0)
; co = internal sound speed: calculated at (dimx/2,dimy/2,dimz/2)
; ce = external sound speed: calculated at (0,0,0)
; aa = radius of cylinder
; r0 = center value of grid along x axis : gridx[-1]/2
; gridx = grid along x axis (radial direction)
; gridy = grid along y axis (radial direction)
; gridz = grid along z axis (longitudinal direction)
; dimx = number of points along x axis
; dimy = number of points along y axis
; dimz = number of points along z axis
; tarr = time array
; wk_rt = w/k solutions corresponding to trapped modes
; ka_rt = ka values of w/k solutions corresponding to trapped modes
; kafix = location in ka_rt array of the ka value considered
; te = (y,z) array of temperature at position ix and time it
; rho = (y,z) array of density at position ix and time it
; vr = (y,z) array of radial velocity at position ix and time it
; vz = (y,z) array of longitudinal velocity at position ix and time it
; velx = (y,z) array of velocity along x at position ix and time it
; vely = (y,z) array of velocity along y at position ix and time it
; velz = (y,z) array of velocity along z at position ix and time it
; NOTE: 
; If 2D model (set = 25) then only quantities with x and z values are returned

;if set eq 2 then begin dir = '/volume1/scratch/set2/' & kanm = 'ka2.24' & endif
;if set eq 3 then begin dir = '/volume1/scratch/set3/' & kanm = 'ka1.27' & endif

if set eq 2 then begin dir = '../../cubes/set2/' & kanm = 'ka2.24' & endif
if set eq 24 then begin dir = '../../cubes/set2/' & kanm = 'ka2.24_hgres' & endif
if set eq 21 then begin dir = '../../cubes/set2/' & kanm = 'ka2.24_highT' & endif
if set eq 22 then begin dir = '../../cubes/set2/' & kanm = 'ka2.24_sm' & endif
if set eq 25 then begin dir = '../../cubes/set2/' & kanm = 'ka2.24_hgres2d' & endif
if set eq 3 then begin dir = '../../cubes/set3/' & kanm = 'ka1.25' & endif
;snum = fix(strmid(slize,11,3))
;if set eq 24 then begin
;   dimx = 2040
;   dimy = 2040
;   dimz = 1142
;   dimt = 30
;   z_u = 27.9875
;   gridx = findgen(dimx)/float(dimx-1)*50.
;   gridy = findgen(dimy)/float(dimy-1)*50.
;   gridz = findgen(dimz)/(dimz-1.)*z_u
;   r0 = gridx[dimx-1]/2.
;   t_u = 4.57641
;   tarr = findgen(dimt)/(dimt-1.)*t_u
;   aa = 10.
;endif

;   restore,'cubes_'+string(it,format='(i3.3)')+'.sav'

if set ne 25 then begin
   restore,dir+'params_'+kanm+'.sav'
   if set ne 24 then begin
      restore, dir+'slice_'+kanm+'_'+string(it,format='(i3.3)')+'t_'+string(ix,format='(i3.3)')+'x'+'.sav'
   endif else begin
      restore, dir+'slice_rh_'+kanm+'_'+string(it,format='(i3.3)')+'t_'+string(ix,format='(i4.4)')+'x'+'.sav'
      restore, dir+'slice_te_'+kanm+'_'+string(it,format='(i3.3)')+'t_'+string(ix,format='(i4.4)')+'x'+'.sav'
   endelse
   gridx0=gridx[ix]
;   gridx = gridx0[ix]
   sizes = size(te)
   dimx0 = sizes[1] & dimy = sizes[2] & dimz = sizes[3]
   dimx = dimy
   if set ne 24 then begin
      velx = fltarr(dimx0,dimy,dimz)
      vely = fltarr(dimx0,dimy,dimz)
      velz = fltarr(dimx0,dimy,dimz)
      for j=0,dimy-1 do begin
         velx[0,j,*] = vr[0,j,*]*(gridx0[0]-r0)/sqrt((gridx0[0]-r0)^2+(gridy[j]-r0)^2)
         vely[0,j,*] = vr[0,j,*]*(gridy[j]-r0)/sqrt((gridx0[0]-r0)^2+(gridy[j]-r0)^2)
      endfor
   endif
endif else begin
   restore,dir+'params_'+kanm+'.sav'
   restore, dir+'cubes_'+kanm+'_'+string(it,format='(i3.3)')+'.sav'
   sizes = size(te_cube)
   dimy = sizes[1] & dimz = sizes[2]
   te = te_cube & rho = rh_cube & vr = vr_cube & vz = vz_cube
   vely = vr*0. & gridy = gridx
   for j=0,dimy-1 do vely[j,*] = vr[j,*]*(gridy[j]-r0)/sqrt((gridy[j]-r0)^2)
   velz = vz_cube
endelse

end
