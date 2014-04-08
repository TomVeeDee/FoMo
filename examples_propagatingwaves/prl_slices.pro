
pro prl_slices,set=set, mua_d=mua_d, ion=ion, w0=w0, imaging=imaging, dx=dx, dy=dy, dz=dz, wayemi=wayemi

; Wrapper routine that produces intensity cubes. Checks for existence
; of previously produced cubes before execution. In this way it can be ran in
; several cores in parallel. 
; INPUT:
; set: = defines model:
;     = 2: 'base model - 171' corresponds to phase speed 1850 km/s, line Fe IX 171
;     = 21: 'high T model' base model with high external temperature
;     = 22: 'smooth model' base model with smooth density profile across cylinder
;     = 23: 'base model - 193' base model with line Fe XII 193
;     = 3 : 'long lambda model' corresponds to ka = 1.25
;     = 24: 'high_res' corresponds to base model with high resolution
;     = 25: 'high_res2d' corresponds to 2d base model with high resolution
; mua_d = float array containing line-of-sight angles: mua_d=[0.,30.,45.,60.]
; ion = string containing the name of the ion: 
;       'fe_9' for Fe IX, 'fe_12' for Fe XII
; w0 = central wavelength of line
; dx = spatial grid size along x: (gridx[dimx-1]-gridx[0])/dimx
; dy = spatial grid size along y: (gridy[dimy-1]-gridy[0])/dimy
; dz = spatial grid size along z: (gridz[dimz-1]-gridz[0])/dimz
 
  if set eq 3 then dir = '/users/cpa/sgijsen/FoMo/examples_propagating/Run_t41_propsauP2L15/'
  if set eq 2 or set eq 21 or set eq 22 or set eq 25 then dir = '/users/cpa/sgijsen/FoMo/examples_propagating/Run_t41_propsauP2L15/'
;  if set eq 24 then dir = '/volume1/scratch/set2/hgres/'
  if set eq 24 then dir = '/users/cpa/pantolin/Modeling/cubes/set2/'

  if set eq 25 then files = file_search(dir+'cubes_*.sav',count=nfiles,/fully_qualify_path) else  files = file_search(dir+'slice_*.sav',count=nfiles,/fully_qualify_path)
  if set eq 24 then begin 
     files0 = file_search(dir+'slice_rh*.sav',count=nfiles0,/fully_qualify_path)
     files1 = file_search(dir+'slice_te*.sav',count=nfiles1,/fully_qualify_path)
  endif
;  files = file_search(dir+'slice_*_014t_*.sav',count=nfiles,/fully_qualify_path)
  if set eq 24 then begin & files = files0 & nfiles = nfiles0 & endif
  restore,files[0]

  params, set=set, it=0, ix=0, ro=ro, re=re, valfven=va, vae=vae, co=co, ce=ce, aa=aa, r0=r0, gridx=gridx, gridy=gridy, gridz=gridz, dimx=dimx, dimy=dimy, dimz=dimz, tarr=tarr, kafix=kafix, ka_rt=ka_rt, wk_rt=wk_rt, te=te, rho=rho, vr=vr, vz=vz, velx=velx, vely=vely, velz=velz

  siz = size(rho)
  nwave = 100
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 or set eq 25 then ndz = 3
  if set eq 3 then ndz = 2
  nsl = 6
  if siz[0] eq 3 then begin
     gridz_ext = fltarr(siz[3]*ndz)
     if imaging eq 0 then begin
        velx_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
        vely_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
        velz_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
        emission_goft_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
     endif else begin 
        intens_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
     endelse
     gridx = gridy
  endif else begin
     gridz_ext = fltarr(siz[2]*ndz)
     if imaging eq 0 then begin
        vely_ext = fltarr(siz[1],siz[2]*ndz)
        velz_ext = fltarr(siz[1],siz[2]*ndz)
        emission_goft_ext = fltarr(siz[1],siz[2]*ndz) 
     endif else begin
        intens_ext = fltarr(siz[1],siz[2]*ndz)
     endelse
  endelse
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 or set eq 25 then begin
     gridz_ext[0:dimz-1] = gridz
     gridz_ext[dimz:2*dimz-1] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
     gridz_ext[dimz*2:*] = gridz_ext[dimz:dimz*2-1]+gridz[dimz-1]+gridz[1]
  endif
  if set eq 3 then begin
     gridz_ext[0:dimz-1] = gridz
     gridz_ext[dimz:*] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
  endif
  if imaging eq 1 then begin
     velz_ext = 0. & vely_ext = 0.
     for j=0,n_elements(mua_d)-1 do begin
        mua = mua_d[j]
        gridlos, gridx=gridz_ext, gridy=gridy, mua_d=mua, velx=velz_ext, vely=vely_ext, dx=dz, dy=dy, n_gridz_ext, n_gridy_ext, ngrid_ext, dl=dl, losvel_ext
        if mua eq 0. then begin 
           n_gridz_ext_00 = n_gridz_ext
           n_gridy_ext_00 = n_gridy_ext
           ngrid_ext_00 = ngrid_ext
        endif
        if mua eq 30. then begin 
           n_gridz_ext_30 = n_gridz_ext
           n_gridy_ext_30 = n_gridy_ext
           ngrid_ext_30 = ngrid_ext
        endif
        if mua eq 45. then begin 
           n_gridz_ext_45 = n_gridz_ext
           n_gridy_ext_45 = n_gridy_ext
           ngrid_ext_45 = ngrid_ext
        endif
        if mua eq 60. then begin 
           n_gridz_ext_60 = n_gridz_ext
           n_gridy_ext_60 = n_gridy_ext
           ngrid_ext_60 = ngrid_ext
        endif        
     endfor
  endif
  for i=0,nfiles-1 do begin
     checkfile=files[i]+'.'+ion+'.done'
     slize = (strsplit(files[i],'/',/extract))[nsl]
     if set eq 2 or set eq 3 then begin
        tmp = strmid(slize,13,3)
        ntsl = fix(tmp)
        tmp = strmid(slize,18,3)
        nxsl = fix(tmp)
     endif
     if set eq 21 then begin
        tmp = strmid(slize,19,3)
        ntsl = fix(tmp)
        tmp = strmid(slize,24,3)
        nxsl = fix(tmp)
     endif
     if set eq 24 then begin
        tmp = strmid(slize,22,3)
        ntsl = fix(tmp)
        tmp = strmid(slize,27,4)
        nxsl = fix(tmp)
     endif
     if set eq 22 then begin
        tmp = strmid(slize,16,3)
        ntsl = fix(tmp)
        tmp = strmid(slize,21,3)
        nxsl = fix(tmp)
     endif
     if set eq 25 then begin
        tmp = strmid(slize,21,3)
        ntsl = fix(tmp)
        nxsl = 0
     endif
     if (file_test(checkfile) eq 0) then begin
        OPENW, lun, checkfile, /get_lun
        free_lun, lun
        print,'Doing slice:',slize
        params, set=set, it=ntsl, ix=nxsl, ro=ro, re=re, valfven=va, vae=vae, co=co, ce=ce, aa=aa, r0=r0, gridx=gridx, gridy=gridy, gridz=gridz, dimx=dimx, dimy=dimy, dimz=dimz, tarr=tarr, kafix=kafix, ka_rt=ka_rt, wk_rt=wk_rt, te=te, rho=rho, vr=vr, vz=vz, velx=velx, vely=vely, velz=velz
        if imaging eq 1 then begin
           lineongrid_int_tab, rho, te, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, logt=logt, watom=watom,intens=intens
        endif else begin
           lineongrid_goft_tab, rho, te, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, emission_goft=emission_goft, goft=goft, logt=logt,wayemi=wayemi, watom=watom, conv=conv,file_abund=file_abund,vers=vers
        endelse
        if set eq 2 or set eq 21 or set eq 22 or set eq 24 then begin
           if imaging eq 0 then begin
              emission_goft_ext[*,*,0:dimz-1] = emission_goft
              emission_goft_ext[*,*,dimz:2*dimz-1] = emission_goft
              emission_goft_ext[*,*,dimz*2:*] = emission_goft
              velx_ext[*,*,0:dimz-1] = velx 
              velx_ext[*,*,dimz:2*dimz-1] = velx 
              velx_ext[*,*,dimz*2:*] = velx 
              vely_ext[*,*,0:dimz-1] = vely
              vely_ext[*,*,dimz:2*dimz-1] = vely
              vely_ext[*,*,dimz*2:*] = vely
              velz_ext[*,*,0:dimz-1] = velz 
              velz_ext[*,*,dimz:2*dimz-1] = velz
              velz_ext[*,*,dimz*2:*] = velz
           endif else begin
              intens_ext[*,*,0:dimz-1] = intens
              intens_ext[*,*,dimz:2*dimz-1] = intens
              intens_ext[*,*,dimz*2:*] = intens
              velx_ext = 0. & vely_ext = 0. & velz_ext = 0.
           endelse
        endif
        if set eq 3 then begin
           if imaging eq 0 then begin
              emission_goft_ext[*,*,0:dimz-1] = emission_goft
              emission_goft_ext[*,*,dimz:*] = emission_goft
              velx_ext[*,*,0:dimz-1] = velx 
              velx_ext[*,*,dimz:*] = velx 
              vely_ext[*,*,0:dimz-1] = vely
              vely_ext[*,*,dimz:*] = vely
              velz_ext[*,*,0:dimz-1] = velz 
              velz_ext[*,*,dimz:*] = velz
           endif else begin
              intens_ext[*,*,0:dimz-1] = intens
              intens_ext[*,*,dimz:*] = intens
              velx_ext = 0. & vely_ext = 0. & velz_ext = 0.
           endelse
        endif
        if set eq 25 then begin
           if imaging eq 0 then begin
              emission_goft_ext[*,0:dimz-1] = emission_goft
              emission_goft_ext[*,dimz:2*dimz-1] = emission_goft
              emission_goft_ext[*,dimz*2:*] = emission_goft
              vely_ext[*,0:dimz-1] = vely
              vely_ext[*,dimz:2*dimz-1] = vely
              vely_ext[*,dimz*2:*] = vely
              velz_ext[*,0:dimz-1] = velz 
              velz_ext[*,dimz:2*dimz-1] = velz
              velz_ext[*,dimz*2:*] = velz
           endif else begin
              intens_ext[*,0:dimz-1] = intens
              intens_ext[*,dimz:2*dimz-1] = intens
              intens_ext[*,dimz*2:*] = intens
              vely_ext = 0. & velz_ext = 0.
           endelse
           gridy = gridx
        endif
        for j=0,n_elements(mua_d)-1 do begin
           mua = mua_d[j]
           if imaging eq 0 then gridlos, gridx=gridz_ext, gridy=gridy, mua_d=mua, velx=velz_ext, vely=vely_ext, dx=dz, dy=dy, n_gridz_ext, n_gridy_ext, ngrid_ext, dl=dl, losvel_ext
           if imaging eq 1 then begin
              if mua eq 0. then begin 
                 n_gridz_ext = n_gridz_ext_00
                 n_gridy_ext = n_gridy_ext_00
                 ngrid_ext = ngrid_ext_00
              endif
              if mua eq 30. then begin 
                 n_gridz_ext = n_gridz_ext_30
                 n_gridy_ext = n_gridy_ext_30
                 ngrid_ext = ngrid_ext_30
              endif
              if mua eq 45. then begin 
                 n_gridz_ext = n_gridz_ext_45
                 n_gridy_ext = n_gridy_ext_45
                 ngrid_ext = ngrid_ext_45
              endif
              if mua eq 60. then begin 
                 n_gridz_ext = n_gridz_ext_60
                 n_gridy_ext = n_gridy_ext_60
                 ngrid_ext = ngrid_ext_60
              endif
           endif
           if mua eq 0. then direction = 2 else direction = 4
           if imaging eq 0 then begin
              integrateemission,emission=emission_goft_ext,logt=logt,n_gridx=n_gridy_ext,n_gridy=n_gridz_ext,ngrid=ngrid_ext,wave=wave,w0=w0,direction=direction,losvel=losvel_ext,imaging=imaging,watom=watom,image,wayemi=wayemi
           endif else begin
              integrateemission,emission=intens_ext,logt=logt,n_gridx=n_gridy_ext,n_gridy=n_gridz_ext,ngrid=ngrid_ext,wave=wave,w0=w0,direction=direction,losvel=losvel_ext,imaging=imaging,image
           endelse
           if mua eq 0.  then image00d_ext=image
           if mua eq 30. then image30d_ext=image
           if mua eq 45. then image45d_ext=image
           if mua eq 60. then image60d_ext=image
        endfor
        if set ne 24 then begin 
           savn0 = (strsplit(slize,'.',/extract)) 
           savn = savn0[0]+'.'+savn0[1]+'.'+ion
        endif else begin
           savn0 = (strsplit(slize,'_',/extract))
           savnp0 = (strsplit(savn0[5],'.',/extract))
           savn = savn0[0]+'_'+savn0[2]+'_'+savn0[3]+'_'+savn0[4]+'_'+savnp0[0]+'.'+ion
        endelse
        if imaging eq 0 then begin
           save,wave,gridx,gridy,gridz_ext,emission_goft,image00d_ext,image30d_ext,image45d_ext,image60d_ext,filename=dir+'rslt_'+savn+'.sav' 
        endif else begin
           save,wave,gridx,gridy,gridz_ext,image00d_ext,image30d_ext,image45d_ext,image60d_ext,filename=dir+'rslt_'+savn+'.sav'
        endelse
     endif
  endfor

end
     
