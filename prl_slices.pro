
pro prl_slices,set=set

; Wrapper routine that produces intensity cubes. Checks for existence
; of previously produced cubes before execution. In this way it can be ran in
; several cores in parallel. 

;  dir = '/volume1/scratch/set2/'
;  set = 24
  if set eq 3 then dir = '/users/cpa/pantolin/Modeling/cubes/set3/'
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 then dir = '/users/cpa/pantolin/Modeling/cubes/set2/'

;  files = file_search(dir+'slice_*.sav',count=nfiles,/fully_qualify_path)
  files = file_search(dir+'slice_*_014t_*.sav',count=nfiles,/fully_qualify_path)
  direction = 4
  mua_dd=[0.,30.]
;  mua_dd=[0.,30.,45.,60.]
  ion = 'fe_9'
;  ion = 'fe_12'

  params, set=set, it=0, ix=0, ro=ro, re=re, va=va, vae=vae, co=co, ce=ce, aa=aa, r0=r0, gridx=gridx, gridy=gridy, gridz=gridz, dimx=dimx, dimy=dimy, dimz=dimz, tarr=tarr, kafix=kafix, ka_rt=ka_rt, wk_rt=wk_rt, te=te, rho=rho, vr=vr, vz=vz, velx=velx, vely=vely, velz=velz

  restore,files[0]
  siz = size(rho)
  nwave = 100
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 then ndz = 3
  if set eq 3 then ndz = 2
  gridz_ext = fltarr(siz[3]*ndz)
  emission_goft_ext = fltarr(siz[1],siz[2],siz[3]*ndz,nwave)
  rho_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
  velx_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
  vely_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
  velz_ext = fltarr(siz[1],siz[2],siz[3]*ndz)
  if set eq 2 or set eq 21 or set eq 22 or set eq 24 then begin
     gridz_ext[0:dimz-1] = gridz
     gridz_ext[dimz:2*dimz-1] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
     gridz_ext[dimz*2:-1] = gridz_ext[dimz:dimz*2-1]+gridz[dimz-1]+gridz[1]
  endif
  if set eq 3 then begin
     gridz_ext[0:dimz-1] = gridz
     gridz_ext[dimz:-1] = gridz[0:dimz-1]+gridz[dimz-1]+gridz[1]
  endif
  gridx = gridy
  for i=0,nfiles-1 do begin
     checkfile=files[i]+'.done'
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
     if (file_test(checkfile) eq 0) then begin
        OPENW, lun, checkfile, /get_lun
        free_lun, lun
        print,'Doing slice:',slize
        params, set, ntsl,nxsl, ro, re, va, vae, co, ce, aa, r0, gridx, gridy, gridz, dimx, dimy, dimz, tarr, kafix, ka_rt, wk_rt, te, rho, vr, vz, velx, vely, velz
        lineongrid, rho, te, wave=wave,nwave=nwave,w0=w0,minwave=minwave,maxwave=maxwave,ion=ion, emission=emission_goft,outwave=outwave,wayemi=2
        if set eq 2 or set eq 21 or set eq 22 or set eq 24 then begin
           emission_goft_ext[*,*,0:dimz-1,*] = emission_goft
           emission_goft_ext[*,*,dimz:2*dimz-1,*] = emission_goft
           emission_goft_ext[*,*,dimz*2:-1,*] = emission_goft
           velx_ext[*,*,0:dimz-1] = velx 
           velx_ext[*,*,dimz:2*dimz-1] = velx 
           velx_ext[*,*,dimz*2:-1] = velx 
           vely_ext[*,*,0:dimz-1] = vely
           vely_ext[*,*,dimz:2*dimz-1] = vely
           vely_ext[*,*,dimz*2:-1] = vely
           velz_ext[*,*,0:dimz-1] = velz 
           velz_ext[*,*,dimz:2*dimz-1] = velz
           velz_ext[*,*,dimz*2:-1] = velz
           rho_ext[*,*,0:dimz-1] = rho 
           rho_ext[*,*,dimz:2*dimz-1] = rho
           rho_ext[*,*,dimz*2:-1] = rho
        endif
        if set eq 3 then begin
           emission_goft_ext[*,*,0:dimz-1,*] = emission_goft
           emission_goft_ext[*,*,dimz:-1,*] = emission_goft
           velx_ext[*,*,0:dimz-1] = velx 
           velx_ext[*,*,dimz:-1] = velx 
           vely_ext[*,*,0:dimz-1] = vely
           vely_ext[*,*,dimz:-1] = vely
           velz_ext[*,*,0:dimz-1] = velz 
           velz_ext[*,*,dimz:-1] = velz
           rho_ext[*,*,0:dimz-1] = rho 
           rho_ext[*,*,dimz:-1] = rho
        endif
        gridx = gridy
        for j=0,n_elements(mua_dd)-1 do begin
           mua_d = mua_dd[j]
           gridlos, gridx=gridz_ext, gridy=gridy, mua_d=mua_d, velx=velz_ext, vely=vely_ext, n_gridz_ext, n_gridy_ext, ngrid_ext, losvel_ext
           integrateemission,emis=emission_goft_ext,n_gridx=n_gridz_ext,n_gridy=n_gridy_ext,ngrid=ngrid_ext,wave=wave,w0=w0,direction=direction,losvel=losvel_ext,image
           if j eq 0 then image00d_ext=image
           if j eq 1 then image30d_ext=image
           if j eq 2 then image45d_ext=image
           if j eq 3 then image60d_ext=image
        endfor
        savn0 = (strsplit(slize,'.',/extract))
        savn = savn0[0]+'.'+savn0[1]+'.'+ion
        save,wave,gridx,gridy,gridz_ext,emission_goft,image00d_ext,image30d_ext,filename=dir+'rslt_'+savn+'.sav'
;        save,wave,gridx,gridy,gridz_ext,emission_goft,image00d_ext,image30d_ext,image45d_ext,image60d_ext,filename=dir+'rslt_'+savn+'.sav'
     endif
  endfor

end
     
