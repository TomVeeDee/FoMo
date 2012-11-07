
pro prl_slices,set=set, mua_d=mua_d, ion=ion, imaging=imaging

; Wrapper routine that produces intensity cubes. Checks for existence
; of previously produced cubes before execution. In this way it can be ran in
; several cores in parallel. 
; INPUT:
; set: = defines model:
;     = 2: 'base model - 171' corresponds to ka = 2.24, line Fe IX 171
;     = 21: 'high T model' base model with high external temperature
;     = 22: 'smooth model' base model with smooth density profile across cylinder
;     = 23: 'base model - 193' base model with line Fe XII 193
;     = 3 : 'long lambda model' corresponds to ka = 1.25
;     = 24: 'high_res' corresponds to base model with high resolution
;     = 25: 'high_res2d' corresponds to 2d base model with high resolution
; mua_d = float array containing line-of-sight angles: mua_d=[0.,30.,45.,60.]
; ion = string containing the name of the ion: 
;       'fe_9' for Fe IX 171, 'fe_12' for Fe XII 193

  if set eq 3 then dir = '/users/cpa/pantolin/Modeling/cubes/set3/'
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 or set eq 25 then dir = '/users/cpa/pantolin/Modeling/cubes/set2/'

  if set eq 25 then files = file_search(dir+'cubes_*.sav',count=nfiles,/fully_qualify_path) else  files = file_search(dir+'slice_*.sav',count=nfiles,/fully_qualify_path)
;  files = file_search(dir+'slice_*_014t_*.sav',count=nfiles,/fully_qualify_path)
  params, set=set, it=0, ix=0, ro=ro, re=re, valfven=va, vae=vae, co=co, ce=ce, aa=aa, r0=r0, gridx=gridx, gridy=gridy, gridz=gridz, dimx=dimx, dimy=dimy, dimz=dimz, tarr=tarr, kafix=kafix, ka_rt=ka_rt, wk_rt=wk_rt, te=te, rho=rho, vr=vr, vz=vz, velx=velx, vely=vely, velz=velz

  restore,files[0]
  siz = size(rho)
  nwave = 100
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 or set eq 25 then ndz = 3
  if set eq 3 then ndz = 2
  if siz[0] eq 3 then begin
     gridz_ext = fltarr(siz[3]*ndz)
     if imaging eq 0 then emission_goft_ext = fltarr(siz[1],siz[2],siz[3]*ndz,nwave) else intens_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
     velx_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
     vely_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
     velz_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
     gridx = gridy
  endif else begin
     gridz_ext = fltarr(siz[2]*ndz)
     if imaging eq 0 then emission_goft_ext = fltarr(siz[1],siz[2]*ndz,nwave) else intens_ext = fltarr(siz[1],siz[2]*ndz)
     vely_ext = fltarr(siz[1],siz[2]*ndz)
     velz_ext = fltarr(siz[1],siz[2]*ndz)
  endelse
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 or set eq 25 then begin
     gridz_ext[0:dimz-1] = gridz
     gridz_ext[dimz:2*dimz-1] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
     gridz_ext[dimz*2:-1] = gridz_ext[dimz:dimz*2-1]+gridz[dimz-1]+gridz[1]
  endif
  if set eq 3 then begin
     gridz_ext[0:dimz-1] = gridz
     gridz_ext[dimz:-1] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
  endif
  for i=0,nfiles-1 do begin
     checkfile=files[i]+'.'+ion+'.done'
     slize = (strsplit(files[i],'/',/extract))[6]
     if set eq 2 or set eq 3 then begin
        tmp = strmid(slize,13,3)
        ntsl = fix(tmp)
        tmp = strmid(slize,18,3)
        nxsl = fix(tmp)
     endif
     if set eq 21 or set eq 24 then begin
        tmp = strmid(slize,19,3)
        ntsl = fix(tmp)
        tmp = strmid(slize,24,3)
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
           lineongrid_int_tab, rho, te, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, intens=intens
        endif else begin
           lineongrid_goft_tab, rho, te, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, emission_goft=emission_goft,wayemi=2
        endelse
        if set eq 2 or set eq 21 or set eq 22 or set eq 24 then begin
           if imaging eq 0 then begin
              emission_goft_ext[*,*,0:dimz-1,*] = emission_goft
              emission_goft_ext[*,*,dimz:2*dimz-1,*] = emission_goft
              emission_goft_ext[*,*,dimz*2:*,*] = emission_goft
           endif else begin
              intens_ext[*,*,0:dimz-1] = intens
              intens_ext[*,*,dimz:2*dimz-1] = intens
              intens_ext[*,*,dimz*2:*] = intens
           endelse
           velx_ext[*,*,0:dimz-1] = velx 
           velx_ext[*,*,dimz:2*dimz-1] = velx 
           velx_ext[*,*,dimz*2:-1] = velx 
           vely_ext[*,*,0:dimz-1] = vely
           vely_ext[*,*,dimz:2*dimz-1] = vely
           vely_ext[*,*,dimz*2:-1] = vely
           velz_ext[*,*,0:dimz-1] = velz 
           velz_ext[*,*,dimz:2*dimz-1] = velz
           velz_ext[*,*,dimz*2:-1] = velz
        endif
        if set eq 3 then begin
           if imaging eq 0 then begin
              emission_goft_ext[*,*,0:dimz-1,*] = emission_goft
              emission_goft_ext[*,*,dimz:*,*] = emission_goft
           endif else begin
              intens_ext[*,*,0:dimz-1,*] = intens
              intens_ext[*,*,dimz:*,*] = intens
           endelse
           velx_ext[*,*,0:dimz-1] = velx 
           velx_ext[*,*,dimz:-1] = velx 
           vely_ext[*,*,0:dimz-1] = vely
           vely_ext[*,*,dimz:-1] = vely
           velz_ext[*,*,0:dimz-1] = velz 
           velz_ext[*,*,dimz:-1] = velz
        endif
        if set eq 25 then begin
           if imaging eq 0 then begin
              emission_goft_ext[*,0:dimz-1,*] = emission_goft
              emission_goft_ext[*,dimz:2*dimz-1,*] = emission_goft
              emission_goft_ext[*,dimz*2:-1,*] = emission_goft
           endif else begin
              intens_ext[*,0:dimz-1,*] = intens
              intens_ext[*,dimz:2*dimz-1,*] = intens
              intens_ext[*,dimz*2:-1,*] = intens
           endelse
           vely_ext[*,0:dimz-1] = vely
           vely_ext[*,dimz:2*dimz-1] = vely
           vely_ext[*,dimz*2:-1] = vely
           velz_ext[*,0:dimz-1] = velz 
           velz_ext[*,dimz:2*dimz-1] = velz
           velz_ext[*,dimz*2:-1] = velz
           gridy = gridx
        endif
        for j=0,n_elements(mua_d)-1 do begin
           mua = mua_d[j]
           gridlos, gridx=gridz_ext, gridy=gridy, mua_d=mua, velx=velz_ext, vely=vely_ext, n_gridz_ext, n_gridy_ext, ngrid_ext, losvel_ext
           if mua eq 0. then direction = 1 else direction = 4
           if imaging eq 0 then begin
              integrateemission,emission=emission_goft_ext,n_gridx=n_gridz_ext,n_gridy=n_gridy_ext,ngrid=ngrid_ext,wave=wave,w0=w0,direction=direction,losvel=losvel_ext,imaging=imaging,image
           endif else begin
              integrateemission,emission=intens,n_gridx=n_gridz_ext,n_gridy=n_gridy_ext,ngrid=ngrid_ext,wave=wave,w0=w0,direction=direction,losvel=losvel_ext,imaging=imaging,image
           endelse
           if mua eq 0.  then image00d_ext=image
           if mua eq 30. then image30d_ext=image
           if mua eq 45. then image45d_ext=image
           if mua eq 60. then image60d_ext=image
        endfor
        savn0 = (strsplit(slize,'.',/extract))
        savn = savn0[0]+'.'+savn0[1]+'.'+ion
        if imaging eq 0 then begin
           save,wave,gridx,gridy,gridz_ext,emission_goft,image00d_ext,image30d_ext,image45d_ext,image60d_ext,filename=dir+'rslt_'+savn+'.sav' 
        endif else begin
           save,wave,gridx,gridy,gridz_ext,intens,image00d_ext,image30d_ext,image45d_ext,image60d_ext,filename=dir+'rslt_'+savn+'.sav'
     endif
  endfor

end
     
