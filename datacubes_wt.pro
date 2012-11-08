
pro datacubes_wt, rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, radius = aa, gridx=gridx, gridz=gridz, dimt=dimt, tarr=tarr, ka_rt=ka_rt, kafix=kafix, wk_rt=wk_rt, vr_t=vr_t, vz_t=vz_t, rtot_t=rtot_t, te_t=te_t, model=model, vr_cube, vz_cube, te_cube, rh_cube, sngcub=sngcub, smooth=smooth, save_cubes=save_cubes

if n_params(0) lt 1 then begin
   print,'Check output directories first'
   print,'datacubes_wt,rho_int=ro, rho_ext=re, valfv_int=va, valv_ext=vae, cs_int=co, cs_ext=ce, radius=aa, gridx=gridx, gridz=gridz, dimt=dimt, tarr=tarr, ka_rt=ka_rt, kafix=kafix, wk_rt=wk_rt, vr_t=vr_t, vz_t=vz_t, rtot_t=rtot_t, te_t=te_t, model=model, vr_cube, vz_cube, te_cube, rh_cube [,sngcub=sngcub,smooth=smooth, save_cubes=save_cubes]'
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
; model = string with name of treated model:
;      model = 'base'corresponds to ka = 2.24, 
;      model = 'long' corresponds to ka = 1.25
;      model = 'high_res' corresponds to ka = 2.24, high spatial resolution
;      model = 'high_res2d' corresponds to ka = 2.24, high 2D spatial resolution
;      model = 'highT' corresponds to ka = 2.24, high external
;      temperature
; OUTPUT:
; Produces IDL save files of cubes for each time step
; vr_cube = vr_cube(dimx,(dimy),dimz) : radial velocity at last time step
; vz_cube = vz_cube(dimx,(dimy),dimz) : longitudinal velocity at last time step
; te_cube = te_cube(dimx,(dimy),dimz) : temperature at last time step
; rh_cube = rh_cube(dimx,(dimy),dimz) : total density at last time step
; OPTIONAL:
; if keyword 'smooth' is set then smooth density profile is used instead
; of sharp profile. In this case run 'smprof_cubes.pro' first and set
; 'rtot_t=rtot_sm_t, te_t=te_sm_t' in call
; if keyword 'save' is set then saves the cubes at each time step
; together with a parameter list file

if keyword_set(smooth) then sm = 1 else sm = 0 
if keyword_set(save_cubes) then sav = 1 else sav = 0
if keyword_set(sngcub) eq 0 then sngcub = 'all'

dimx = n_elements(gridx)
dimz = n_elements(gridz)
r0 = gridx[dimx-1]/2.

; OUTPUT DIRECTORY:
cubedir='../cubes/set2/'

if model eq 'base' then  kanm = 'ka2.24'
if model eq 'long' then  kanm = 'ka1.25'
if model eq 'high_res' then  kanm = 'ka2.24_hgres'
if model eq 'high_res2d' then  kanm = 'ka2.24_hgres2d'
if model eq 'highT' then kanm = 'ka2.24highT'

fixk = 1
if model eq 'ka2.24_hgres2d' then begin
   vr_cube=fltarr(dimx,dimz)
   vz_cube=fltarr(dimx,dimz)
   te_cube=fltarr(dimx,dimz)
   rh_cube=fltarr(dimx,dimz)

   for j=0,dimt-1 do begin
      if sav eq 1 and sm eq 0 then begin
         vr_cube = reform(vr_t[*,*,j])
         vz_cube = reform(vz_t[*,*,j])
         rh_cube = reform(rtot_t[*,*,j])
         te_cube = reform(te_t[*,*,j])
         save,vr_cube,vz_cube,te_cube,rh_cube,filename=cubedir+'cubes2d_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
      endif
      if sav eq 1 and sm eq 1 then begin
         vr_cube_sm = smooth(reform(vr_t[*,*,j]),[4,4])
         vz_cube_sm = smooth(reform(vz_t[*,*,j]),[4,4])
         rh_cube_sm = smooth(reform(rtot_t[*,*,j]),[4,4])
         te_cube_sm = smooth(reform(te_t[*,*,j]),[4,4])
         save,vr_cube_sm,vz_cube_sm,te_cube_sm,rh_cube_sm,filename=cubedir+'cubes_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
      endif
      print,string(13b)+' % finished: ',float(j)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
   if sav eq 1 then save, ro, re, va, vae, co, ce, aa, r0, gridx, gridz, dimx, dimz, dimt, tarr, kafix, ka_rt, wk_rt,filename=cubedir+'params_'+kanm+'.sav'
endif else begin
   gridy = gridx
   dimy = dimx
   print,'doing '+sngcub+' cube'
   if sngcub eq 'all' or sngcub eq 'vr' then vr_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'vz' then vz_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'te' then te_cube=fltarr(dimx,dimy,dimz)
   if sngcub eq 'all' or sngcub eq 'rh' then rh_cube=fltarr(dimx,dimy,dimz)
   dist_cube=fltarr(dimx,dimy,dimz)
   distance = fltarr(dimx,dimy)
   print,'1st step check'
; (frequency, x, y, z, time)
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
            lgrid = ([min(abs(locdis[l]-abs(gridx[0:dimx/2]-r0))),!c])[1]
            if sngcub eq 'all' or sngcub eq 'vr' then col3dvr = fltarr(ndist*dimz)
            if sngcub eq 'all' or sngcub eq 'vz' then col3dvz = fltarr(ndist*dimz)
            if sngcub eq 'all' or sngcub eq 'te' then col3dte = fltarr(ndist*dimz)
            if sngcub eq 'all' or sngcub eq 'rh' then col3drh = fltarr(ndist*dimz)
            if fixk eq 1 then begin
               if sngcub eq 'all' or sngcub eq 'vr' then colvr = reform(vr_t[lgrid,*,j]) 
               if sngcub eq 'all' or sngcub eq 'vz' then colvz = reform(vz_t[lgrid,*,j])
               if sngcub eq 'all' or sngcub eq 'te' then colte = reform(te_t[lgrid,*,j])
               if sngcub eq 'all' or sngcub eq 'rh' then colrh = reform(rtot_t[lgrid,*,j])
            endif else begin
               if sngcub eq 'all' or sngcub eq 'vr' then colvr = reform(vr_t[kafix,lgrid,*,j]) 
               if sngcub eq 'all' or sngcub eq 'vz' then colvz = reform(vz_t[kafix,lgrid,*,j])
               if sngcub eq 'all' or sngcub eq 'te' then colte = reform(te_t[kafix,lgrid,*,j])
               if sngcub eq 'all' or sngcub eq 'rh' then colrh = reform(rtot_t[kafix,lgrid,*,j])
            endelse
            for k=0,dimz-1 do begin
               if sngcub eq 'all' or sngcub eq 'vr' then col3dvr[k*ndist:(k+1)*ndist-1]=replicate(colvr[k],ndist)
               if sngcub eq 'all' or sngcub eq 'vz' then col3dvz[k*ndist:(k+1)*ndist-1]=replicate(colvz[k],ndist)
               if sngcub eq 'all' or sngcub eq 'te' then col3dte[k*ndist:(k+1)*ndist-1]=replicate(colte[k],ndist)
               if sngcub eq 'all' or sngcub eq 'rh' then col3drh[k*ndist:(k+1)*ndist-1]=replicate(colrh[k],ndist)
            endfor
            if sngcub eq 'all' or sngcub eq 'vr' then vr_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dvr
            if sngcub eq 'all' or sngcub eq 'vz' then vz_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dvz
            if sngcub eq 'all' or sngcub eq 'te' then te_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3dte
            if sngcub eq 'all' or sngcub eq 'rh' then rh_cube[Rl[Rl[l]:Rl[l+1]-1]] = col3drh
         endif
      endfor
      if sm eq 1 then begin
         if sngcub eq 'all' or sngcub eq 'vr' then vr_cube_sm = vr_cube 
         if sngcub eq 'all' or sngcub eq 'vz' then vz_cube_sm = vz_cube 
         if sngcub eq 'all' or sngcub eq 'te' then te_cube_sm = te_cube 
         if sngcub eq 'all' or sngcub eq 'rh' then rh_cube_sm = rh_cube 
         for i=0,dimx-1 do begin
            if sngcub eq 'all' or sngcub eq 'vr' then vr_cube_sm[i,*,*]=smooth(reform(vr_cube[i,*,*]),[4,4])
            if sngcub eq 'all' or sngcub eq 'vz' then vz_cube_sm[i,*,*]=smooth(reform(vz_cube[i,*,*]),[4,4])
            if sngcub eq 'all' or sngcub eq 'rh' then rh_cube_sm[i,*,*]=smooth(reform(rh_cube[i,*,*]),[4,4])
            if sngcub eq 'all' or sngcub eq 'te' then te_cube_sm[i,*,*]=smooth(reform(te_cube[i,*,*]),[4,4])
         endfor
      endif
      if sav eq 1 and sm eq 0 then begin
         save,vr_cube,vz_cube,te_cube,rh_cube,filename=cubedir+'cubes_'+sngcub+'_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
      endif
      if sav eq 1 and sm eq 1 then begin
         save,vr_cube_sm,vz_cube_sm,te_cube_sm,rh_cube_sm,filename=cubedir+'cubes_'+sngcub+'_'+kanm+'_'+string(j,format="(i3.3)")+'.sav'
      endif
      print,string(13b)+' % finished: ',float(j)*100./(dimt-1),format='(a,f4.0,$)'
   endfor
   if sav eq 1 then save, ro, re, va, vae, co, ce, aa, r0, gridx, gridy, gridz, dimx, dimy, dimz, dimt, tarr, kafix, ka_rt, wk_rt,filename=cubedir+'params_'+kanm+'.sav'
endelse

end


