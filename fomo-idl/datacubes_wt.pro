
pro datacubes_wt, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, b_int=bo, b_ext=be, radius = aa, gridx=gridx, gridz=gridz, gridr=gridr, dimt=dimt, diml=diml, tarr=tarr, ka_rt=ka_rt, kafix=kafix, wk_rt=wk_rt, vr_t=vr_t, vt_t=vt_t, vz_t=vz_t, rtot_t=rtot_t, te_t=te_t, br_t=br_t, bt_t=bt_t, bz_t=bz_t, btot_t=btot_t, model=model, sngcub=sngcub, mag=mag, save_cubes=save_cubes, nmode=nmode, smooth=smooth

if n_params(0) lt 1 then begin
   print,'Check output directories first'
   print,'datacubes_wt, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, b_int=bo, b_ext=be, radius = aa, gridx=gridx, gridz=gridz, gridr=gridr, dimt=dimt, diml=diml, tarr=tarr, ka_rt=ka_rt, kafix=kafix, wk_rt=wk_rt, vr_t=vr_t, vt_t=vt_t, vz_t=vz_t, rtot_t=rtot_t, te_t=te_t, br_t=br_t, bt_t=bt_t, bz_t=bz_t, btot_t=btot_t, model=model, sngcub=sngcub, mag=mag, save_cubes=save_cubes, nmode=nmode, smooth=smooth'
   return
endif

; Produces 2D or 3D cubes (dimx,(dimy),dimz) of thermodynamic quantities for each
; time step by rotating around the axis of the cylinder
; INPUT:
; ro = internal density: rho[dimx/2,dimy/2,dimz/2]
; re = external density: rho[0,0,0]
; va = internal Alfven velocity: calculated at (dimx/2,dimy/2,dimz/2)
; vae = external Alfven velocity: calculated at (0,0,0)
; co = internal sound speed: calculated at (dimx/2,dimy/2,dimz/2)
; ce = external sound speed: calculated at (0,0,0)
; aa = radius of cylinder
; gridx = grid along x axis (radial direction)
; gridz = grid along z axis (longitudinal direction)
; dimt = number of points in time dimension
; tarr = time array
; wk_rt = w/k solutions corresponding to trapped modes
; ka_rt = ka values of w/k solutions corresponding to trapped modes
; vr_t = vr(dimx,dimz,dimt) : radial velocity
; vz_t = vz(dimx,dimz,dimt) : longitudinal velocity
; rtot_t = rtot(dimx,dimz,dimt) : total density
; te_t = te_t(dimx,dimz,dimt) : temperature
; br_t = br_t(dimx,dimy,dimz) : radial magnetic field
; bz_t = bz_t(dimx,dimy,dimz) : longitudinal magnetic field
; model = string with name of treated model:
;      model = 'base'corresponds to ka = 2.24, 
;      model = 'long' corresponds to ka = 1.25
;      model = 'high_res' corresponds to ka = 2.24, high spatial resolution
;      model = 'high_res2d' corresponds to ka = 2.24, high 2D spatial resolution
;      model = 'highT' corresponds to ka = 2.24, high external
;      temperature
; OUTPUT:
; Produces IDL save files of cubes for each time step
; OPTIONAL:
; if mag keyword is set then produces magnetic field:
; br_cube = te_cube(dimx,(dimy),dimz) : temperature at last time step
; bz_cube = rh_cube(dimx,(dimy),dimz) : total density at last time step
; if keyword 'smooth' is set then smooth density profile is used instead
; of sharp profile. In this case run 'smprof_cubes.pro' first
; if keyword 'save' is set then saves the cubes at each time step
; together with a parameter list file

if keyword_set(smooth) then sm = 1 else sm = 0 
if keyword_set(save_cubes) then sav = 1 else sav = 0
if keyword_set(sngcub) eq 0 then sngcub = 'all'
if keyword_set(nmode) eq 0 then n = 0 else n = nmode

if n eq 0 then mode = 'sausage'
if n eq 1 then mode = 'kink'

dimx = n_elements(gridx)
dimz = n_elements(gridz)
dimr = n_elements(gridr)
r0 = gridx[dimx-1]/2.

; OUTPUT DIRECTORY:
cubedir='../cubes/set2/'

if model eq 'base' then  kanm = 'ka2.24'
if model eq 'long' then  kanm = 'ka1.25'
if model eq 'high_res' then  kanm = 'ka2.24_hgres'
if model eq 'high_res2d' then  kanm = 'ka2.24_hgres2d'
if model eq 'highT' then kanm = 'ka2.24highT'
if model eq 'kink' then kanm = 'kk_ka0.03'

if model eq 'ka2.24_hgres2d' then begin
   vr_cube=fltarr(dimx,dimz)
   vz_cube=fltarr(dimx,dimz)
   te_cube=fltarr(dimx,dimz)
   rh_cube=fltarr(dimx,dimz)
   if keyword_set(mag) then begin
      br_cube=fltarr(dimx,dimz)
      bz_cube=fltarr(dimx,dimz)
   endif

   for j=0,dimt-1 do begin
      if sav eq 1 and sm eq 0 then begin
         vr_cube[0:dimx/2-1,*] = reform(reverse(vr_t[*,*,j],1))
         vr_cube[dimx/2:dimx-1,*] = reform(vr_t[*,*,j])
         vz_cube[0:dimx/2-1,*] = reform(reverse(vz_t[*,*,j],1))
         vz_cube[dimx/2:dimx-1,*] = reform(vz_t[*,*,j])
         rh_cube[0:dimx/2-1,*] = reform(reverse(rh_t[*,*,j],1))
         rh_cube[dimx/2:dimx-1,*] = reform(rh_t[*,*,j])
         te_cube[0:dimx/2-1,*] = reform(reverse(te_t[*,*,j],1))
         te_cube[dimx/2:dimx-1,*] = reform(te_t[*,*,j])
         if keyword_set(mag) then begin
            br_cube[0:dimx/2-1,*] = reform(reverse(br_t[*,*,j],1))
            br_cube[dimx/2:dimx-1,*] = reform(br_t[*,*,j])
            bz_cube[0:dimx/2-1,*] = reform(reverse(bz_t[*,*,j],1))
            bz_cube[dimx/2:dimx-1,*] = reform(bz_t[*,*,j])
         endif
         save,vr_cube,vz_cube,te_cube,rh_cube,br_cube,bz_cube,filename=cubedir+'cubes2d_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
      endif
      if sav eq 1 and sm eq 1 then begin
         vr_cube[0:dimx/2-1,*] = smooth(reform(reverse(vr_t[*,*,j],1)),[4,4])
         vr_cube[dimx/2:dimx-1,*] = smooth(reform(vr_t[*,*,j]),[4,4])
         vz_cube[0:dimx/2-1,*] = smooth(reform(reverse(vz_t[*,*,j],1)),[4,4])
         vz_cube[dimx/2:dimx-1,*] = smooth(reform(vz_t[*,*,j]),[4,4])
         rh_cube[0:dimx/2-1,*] = smooth(reform(reverse(rh_t[*,*,j],1)),[4,4])
         rh_cube[dimx/2:dimx-1,*] = smooth(reform(rh_t[*,*,j]),[4,4])
         te_cube[0:dimx/2-1,*] = smooth(reform(reverse(te_t[*,*,j],1)),[4,4])
         te_cube[dimx/2:dimx-1,*] = smooth(reform(te_t[*,*,j]),[4,4])
         if keyword_set(mag) eq 1 then begin
            br_cube[0:dimx/2-1,*] = smooth(reform(reverse(br_t[*,*,j],1)),[4,4])
            br_cube[dimx/2:dimx-1,*] = smooth(reform(br_t[*,*,j]),[4,4])
            bz_cube[0:dimx/2-1,*] = smooth(reform(reverse(bz_t[*,*,j],1)),[4,4])
            bz_cube[dimx/2:dimx-1,*] = smooth(reform(bz_t[*,*,j]),[4,4])
         endif
         save,vr_cube,vz_cube,te_cube,rh_cube,br_cube,bz_cube,filename=cubedir+'cubes_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
      endif
      print,string(13b)+' % finished: ',float(j)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
   if sav eq 1 then save, ro, re, va, vae, co, ce, bo, be, aa, r0, gridx, gridz, dimx, dimz, dimt, tarr, kafix, ka_rt, wk_rt,filename=cubedir+'params_'+kanm+'.sav'
endif else begin
   gridy = gridx
   dimy = dimx
   print,'doing '+sngcub+' cube'
   if sngcub eq 'all' or sngcub eq 'vr' then vr_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'vt' then vt_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'vz' then vz_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'te' then te_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'rh' then rh_cube=fltarr(dimx,dimy,dimz)
   if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'br')) then br_cube=fltarr(dimx,dimy,dimz)
   if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bt')) then bt_cube=fltarr(dimx,dimy,dimz)
   if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bz')) then bz_cube=fltarr(dimx,dimy,dimz)
   if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'btot')) then btot_cube=fltarr(dimx,dimy,dimz)

   print,'1st step check'
; (frequency, x, y, z, time)
   if mode eq 'sausage' then begin

      dist_cube=fltarr(dimx,dimy,dimz)
      distance = fltarr(dimx,dimy)
      for i=0,dimx-1 do for j=0,dimy-1 do distance[i,j] = sqrt((r0-gridx[i])^2+(r0-gridy[j])^2)
      for i=0,dimz-1 do dist_cube[*,*,i] = distance

      mxdist = sqrt((gridx[dimx-1]-r0)^2+(gridy[dimy-1]-r0)^2)
      if model eq 'ka2.24_hgres' then bsiz = 0.001 else bsiz = 0.01
      cyl = histogram(dist_cube,binsize=bsiz,locations=loccyl,reverse_indices=Rl)
      dis = histogram(dist_cube[*,*,0],binsize=bsiz,locations=locdis)
      ncyl = n_elements(cyl)
      print,'n_elements(cyl):',ncyl
      for j=0,dimt-1 do begin
         for l=0.,ncyl-1 do begin
            ndist = dis[l]
            if ndist gt 0 then begin
;               lgrid = ([min(abs(locdis[l]-abs(gridx[0:dimx/2]-r0))),!c])[1]
               lr = locdis[l]*(dimr-1)/gridr[-1]
               if sngcub eq 'all' or sngcub eq 'vr' then begin colvr = fltarr(dimz) & col3dvr = fltarr(ndist*dimz) & endif
               if sngcub eq 'all' or sngcub eq 'vz' then begin colvz = fltarr(dimz) & col3dvz = fltarr(ndist*dimz) & endif
               if sngcub eq 'all' or sngcub eq 'te' then begin colte = fltarr(dimz) & col3dte = fltarr(ndist*dimz) & endif
               if sngcub eq 'all' or sngcub eq 'rh' then begin colrh = fltarr(dimz) & col3drh = fltarr(ndist*dimz) & endif
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'br')) then begin colbr = fltarr(dimz) & col3dbr = fltarr(ndist*dimz) & endif
;               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bt')) then begin colbt = fltarr(dimz) & col3dbt = fltarr(ndist*dimz) & endif
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bz')) then begin colbz = fltarr(dimz) & col3dbz = fltarr(ndist*dimz) & endif
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'btot')) then begin colbtot = fltarr(dimz) & col3dbtot = fltarr(ndist*dimz) & endif
               if sngcub eq 'all' then begin
                  for ii=0,dimz-1 do begin
                     colvr[ii] = interpolate(reform(vr_t[*,ii,j]),lr)
                     colvz[ii] = interpolate(reform(vz_t[*,ii,j]),lr)
                     colte[ii] = interpolate(reform(te_t[*,ii,j]),lr)
                     colrh[ii] = interpolate(reform(rtot_t[*,ii,j]),lr)
                     if keyword_set(mag) then begin
                        colbr[ii] = interpolate(reform(br_t[*,ii,j]),lr)
;                        colbt[ii] = interpolate(reform(bt_t[*,ii,j]),lr)
                        colbz[ii] = interpolate(reform(bz_t[*,ii,j]),lr)
                        colbtot[ii] = interpolate(reform(btot_t[*,ii,j]),lr)
                     endif
                  endfor
               endif else begin
                  if sngcub eq 'vr' then for ii=0,dimz-1 do colvr[ii] = interpolate(reform(vr_t[*,ii,j]),lr)
                  if sngcub eq 'vz' then for ii=0,dimz-1 do colvz[ii] = interpolate(reform(vz_t[*,ii,j]),lr)
                  if sngcub eq 'te' then for ii=0,dimz-1 do colte[ii] = interpolate(reform(te_t[*,ii,j]),lr)
                  if sngcub eq 'rh' then for ii=0,dimz-1 do colrh[ii] = interpolate(reform(rtot_t[*,ii,j]),lr)
                  if (keyword_set(mag) and sngcub eq 'br') then for ii=0,dimz-1 do colbr[ii] = interpolate(reform(br_t[*,ii,j]),lr)
;                  if (keyword_set(mag) and sngcub eq 'bt') then for ii=0,dimz-1 do colbt[ii] = interpolate(reform(bt_t[*,ii,j]),lr)
                  if (keyword_set(mag) and sngcub eq 'bz') then for ii=0,dimz-1 do colbz[ii] = interpolate(reform(bz_t[*,ii,j]),lr)
                  if (keyword_set(mag) and sngcub eq 'btot') then for ii=0,dimz-1 do colbtot[ii] = interpolate(reform(btot_t[*,ii,j]),lr)
               endelse
               for k=0,dimz-1 do begin
                  if sngcub eq 'all' then begin 
                     col3dvr[k*ndist:(k+1)*ndist-1]=replicate(colvr[k],ndist)
                     col3dvz[k*ndist:(k+1)*ndist-1]=replicate(colvz[k],ndist)
                     col3dte[k*ndist:(k+1)*ndist-1]=replicate(colte[k],ndist)
                     col3drh[k*ndist:(k+1)*ndist-1]=replicate(colrh[k],ndist)
                     if keyword_set(mag) then begin
                        col3dbr[k*ndist:(k+1)*ndist-1]=replicate(colbr[k],ndist)
;                        col3dbt[k*ndist:(k+1)*ndist-1]=replicate(colbt[k],ndist)
                        col3dbz[k*ndist:(k+1)*ndist-1]=replicate(colbz[k],ndist)
                        col3dbtot[k*ndist:(k+1)*ndist-1]=replicate(colbtot[k],ndist)
                     endif
                  endif else begin
                     if sngcub eq 'vr' then col3dvr[k*ndist:(k+1)*ndist-1]=replicate(colvr[k],ndist)
                     if sngcub eq 'vz' then col3dvz[k*ndist:(k+1)*ndist-1]=replicate(colvz[k],ndist)
                     if sngcub eq 'te' then col3dte[k*ndist:(k+1)*ndist-1]=replicate(colte[k],ndist)
                     if sngcub eq 'rh' then col3drh[k*ndist:(k+1)*ndist-1]=replicate(colrh[k],ndist)
                     if (keyword_set(mag) and sngcub eq 'br') then col3dbr[k*ndist:(k+1)*ndist-1]=replicate(colbr[k],ndist)
;                     if (keyword_set(mag) and sngcub eq 'bt') then col3dbt[k*ndist:(k+1)*ndist-1]=replicate(colbt[k],ndist)
                     if (keyword_set(mag) and sngcub eq 'bz') then col3dbz[k*ndist:(k+1)*ndist-1]=replicate(colbz[k],ndist)
                     if (keyword_set(mag) and sngcub eq 'btot') then col3dbtot[k*ndist:(k+1)*ndist-1]=replicate(colbtot[k],ndist)
                  endelse
               endfor
               if sngcub eq 'all' or sngcub eq 'vr' then vr_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dvr
               if sngcub eq 'all' or sngcub eq 'vz' then vz_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dvz
               if sngcub eq 'all' or sngcub eq 'te' then te_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dte
               if sngcub eq 'all' or sngcub eq 'rh' then rh_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3drh
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'br')) then br_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dbr
;               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bt')) then bt_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dbt
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bz')) then bz_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dbz
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'btot')) then btot_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dbtot
            endif
         endfor
         if sm eq 1 then begin
            for i=0,dimx-1 do begin
               if sngcub eq 'all' or sngcub eq 'vr' then vr_cube[i,*,*]=smooth(reform(vr_cube[i,*,*]),[4,4])
               if sngcub eq 'all' or sngcub eq 'vz' then vz_cube[i,*,*]=smooth(reform(vz_cube[i,*,*]),[4,4])
               if sngcub eq 'all' or sngcub eq 'rh' then rh_cube[i,*,*]=smooth(reform(rh_cube[i,*,*]),[4,4])
               if sngcub eq 'all' or sngcub eq 'te' then te_cube[i,*,*]=smooth(reform(te_cube[i,*,*]),[4,4])
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'br')) then br_cube[i,*,*]=smooth(reform(br_cube[i,*,*]),[4,4])
;               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bt')) then bt_cube[i,*,*]=smooth(reform(bt_cube[i,*,*]),[4,4])
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bz')) then bz_cube[i,*,*]=smooth(reform(bz_cube[i,*,*]),[4,4])
               if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'btot')) then btot_cube[i,*,*]=smooth(reform(btot_cube[i,*,*]),[4,4])
            endfor
         endif
         if sav eq 1 then begin
            if keyword_set(mag) then save,vr_cube,vz_cube,te_cube,rh_cube,br_cube,bz_cube,btot_cube,filename=cubedir+'cubes_'+sngcub+'_'+kanm+'_'+string(j,format="(i3.3)")+'.sav' else save,vr_cube,vz_cube,te_cube,rh_cube,filename=cubedir+'cubes_'+sngcub+'_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
         endif
         print,string(13b)+' % finished: ',float(j)*100./(dimt-1),format='(a,f4.0,$)'
      endfor
      if sav eq 1 then save, ro, re, va, vae, co, ce, bo, be, aa, r0, gridx, gridy, gridz, dimx, dimy, dimz, dimt, tarr, kafix, ka_rt, wk_rt,filename=cubedir+'params_'+kanm+'.sav'
   endif
   if mode eq 'kink' then begin
      for k=0,dimt-1 do begin
         for i=0,dimx-1 do begin
            for j=0,dimy-1 do begin
               for l=0,dimz-1 do begin
;                  if i lt dimx/2 then rsig = 1 else rsig = -1
                  r = sqrt((gridx[i]-r0)^2.+(gridy[j]-r0)^2.)
;                  th0 = atan((gridy[j]-r0)/(gridx[i]-r0))+!pi/2.
                  th0 = atan((gridx[j]-r0)/(gridy[i]-r0))+!pi/2.
                  if (gridx[i]-r0) ge 0. and (gridy[j]-r0) ge 0. then th = th0
                  if (gridx[i]-r0) lt 0. and (gridy[j]-r0) gt 0. then th = th0 + !pi
                  if (gridx[i]-r0) lt 0. and (gridy[j]-r0) lt 0. then th = th0 + !pi
                  if (gridx[i]-r0) gt 0. and (gridy[j]-r0) lt 0. then th = th0
;                  lr = (r+r0)*(dimr-1)/(gridr[-1]+r0)
                  lr = r*(dimr-1)/gridr[-1]
                  lth = th*diml/(2*!pi)
                  vr_cube[i,j,l] = interpolate(reform(vr_t[*,*,l,k]),lr,lth)
                  vt_cube[i,j,l] = interpolate(reform(vt_t[*,*,l,k]),lr,lth)
                  vz_cube[i,j,l] = interpolate(reform(vz_t[*,*,l,k]),lr,lth)
                  rh_cube[i,j,l] = interpolate(reform(rtot_t[*,*,l,k]),lr,lth)
                  te_cube[i,j,l] = interpolate(reform(te_t[*,*,l,k]),lr,lth)
                  if (keyword_set(mag)) then begin
                     br_cube[i,j,l] = interpolate(reform(br_t[*,*,l,k]),r,lth)
                     bt_cube[i,j,l] = interpolate(reform(bt_t[*,*,l,k]),r,lth)
                     bz_cube[i,j,l] = interpolate(reform(bz_t[*,*,l,k]),r,lth)
                     btot_cube[i,j,l] = interpolate(reform(btot_t[*,*,l,k]),r,lth)
                  endif
               endfor
            endfor
         endfor
         if sav eq 1 then begin
            if keyword_set(mag) then save,vr_cube,vt_cube,vz_cube,te_cube,rh_cube,br_cube,bt_cube,bz_cube,btot_cube,filename=cubedir+'cubes_'+kanm+'_'+string(k,format="(i3.3)")+'.sav' else save,vr_cube,vt_cube,vz_cube,te_cube,rh_cube,filename=cubedir+'cubes_'+kanm+'_'+string(k,format="(i3.3)")+'.sav'
         endif
         print,string(13b)+' % finished: ',float(k)*100./(dimt-1),format='(a,f4.0,$)'
      endfor
      if sav eq 1 then save, ro, re, va, vae, co, ce, bo, be, aa, r0, gridx, gridy, gridz, dimx, dimy, dimz, dimt, tarr, kafix, ka_rt, wk_rt,filename=cubedir+'params_'+kanm+'.sav'
   endif
endelse

end


