
pro synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,gridz=gridz,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,line_1=line_1,img_1=img_1,line_2=line_2,img_2=img_2,line_3=line_3,img_3=img_3,line_4=line_4,img_4=img_4,conv=conv, wayemi=wayemi,imgfront=imgfront,filenm=filenm,gotdir=gotdir,file_abund=file_abund,vers=vers

if arg_present(rho) lt 1 then begin
   print,'synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,gridz=gridz,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,line_1=line_1,img_1=img_1,line_2=line_2,img_2=img_2,line_3=line_3,img_3=img_3,line_4=line_4,img_4=img_4,conv=conv, wayemi=wayemi,imgfront=imgfront,filenm=filenm,gotdir=gotdir,file_abund=file_abund,vers=vers'
   return
endif

; INPUT:
; rho, tem, v1m, v2m, gridx, gridy, ion, mua_d, w0


if keyword_set(conv) then begin
   proton = 1.67262158*10^(-27.)
   kboltz = 1.380658*10^(-23.)
   normro = 1.e10
   normte = proton/(2*kboltz)*normro
;   nefrac = 1.2
   nefrac = 1.
   rh = rho*normro*nefrac
   te = tem /normte
endif else begin
   ; check for CGS
   rh = rho
   te = tem
   ne_s = nem
endelse

;ion = 'fe_9'
;mua_d = [0., 45., 90.]
direction = 4
imaging = 0
nang = n_elements(mua_d)
width = 0.
;n_e=rh/proton/1.e6 

velx = v1m
vely = v2m
ngridsx = n_elements(gridx)
ngridsy = n_elements(gridy)
dx = (gridx[-1]-gridx[0])/ngridsx
dy = (gridy[-1]-gridy[0])/ngridsy

if wayemi eq 5 then begin
   ; SDO AIA filters:
   imaging = 1
   n_e = nem/1.e6
   logt = alog10(tem)
   if width gt 0. then begin
      dimz = width/dx
      gridz = findgen(dimz)*dx
   endif else begin
      dimz =1
   endelse
   interpol_emiss_data,n_e,tem,ion=ion, w0=w0,emission_goft=emission_goft,filenm=filenm
   imgfront = emission_goft*dimz
   for i=0,nang-1 do begin
      mua = mua_d[i]
      gridlos, gridx=gridx, gridy=gridy, mua_d=mua, velx=velx, vely=vely, dx=dx, dy=dy, n_gridx, n_gridy, ngrid, dl=dl, losvel
      losvel = -losvel/1.e2
      integrateemission,emission=emission_goft,logt=logt,n_gridx=n_gridx,n_gridy=n_gridy,ngrid=ngrid,w0=w0,direction=direction,losvel=losvel,imaging=imaging,imsp
      dlos = (size(imsp))[1]
      inan = string(i,format="(i1)")
      exe1 = 'line_'+inan+' = imsp'
      exe2 = 'dl_'+inan+'= dl'
      exe3 = 'n_gridx_'+inan+'=n_gridx'
      exe4 = 'n_gridy_'+inan+'=n_gridy'
      exe5 = 'ngrid_'+inan+'=ngrid'
      exe6 = 'img_'+inan+'= fltarr(dlos,dimz)'
      exe7 = 'for j=0,dimz-1 do img_'+inan+'[*,j] = line_'+inan

      void = execute(exe1)
      void = execute(exe2)
      void = execute(exe3)
      void = execute(exe4)
      void = execute(exe5)
      if width gt 0. then begin
         void = execute(exe6)
         void = execute(exe7)
      endif
   endfor
endif else begin

   lineongrid_goft_tab, rh, te, ne_s=ne_s, gotdir=gotdir,wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, emission_goft=emission_goft,goft=goft, logt=logt,wayemi=wayemi,watom=watom,conv=conv,file_abund=file_abund,vers=vers

; Frontview:
   if width gt 0. then begin
      dimz = width/dx
      gridz = findgen(dimz)*dx
   endif else begin
      dimz = 1
   endelse
   imgfront = emission_goft*dimz

; Sideview:
   for i=0,nang-1 do begin
      mua = mua_d[i]
      gridlos, gridx=gridx, gridy=gridy, mua_d=mua, velx=velx, vely=vely, dx=dx, dy=dy, n_gridx, n_gridy, ngrid, dl=dl, losvel
      losvel = -losvel/1.e2
      integrateemission,emission=emission_goft,logt=logt,n_gridx=n_gridx,n_gridy=n_gridy,ngrid=ngrid,wave=wave,w0=w0,direction=direction,losvel=losvel,imaging=imaging,imsp,watom=watom,wayemi=wayemi
      
      dlos = (size(imsp))[1]
      inan = string(i,format="(i1)")
      exe1 = 'line_'+inan+' = imsp'
      exe2 = 'dl_'+inan+'= dl'
      exe3 = 'n_gridx_'+inan+'=n_gridx'
      exe4 = 'n_gridy_'+inan+'=n_gridy'
      exe5 = 'ngrid_'+inan+'=ngrid'
      exe6 = 'img_'+inan+'= fltarr(dlos,dimz)'
      exe7 = 'for j=0,dimz-1 do img_'+inan+'[*,j] = line_'+inan

      void = execute(exe1)
      void = execute(exe2)
      void = execute(exe3)
      void = execute(exe4)
      void = execute(exe5)
      if width gt 0. then begin
         void = execute(exe6)
         void = execute(exe7)
      endif
   endfor
endelse

end
