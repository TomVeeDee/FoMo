
pro synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,gridz=gridz,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,line_1=line_1,img_1=img_1,line_2=line_2,img_2=img_2,line_3=line_3,img_3=img_3,line_4=line_4,img_4=img_4,conv=conv, wayemi=wayemi,imgfront=imgfront,filenm=filenm,gotdir=gotdir,file_abund=file_abund,vers=vers

if arg_present(rho) lt 1 then begin
   print,'synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,gridz=gridz,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,line_1=line_1,img_1=img_1,line_2=line_2,img_2=img_2,line_3=line_3,img_3=img_3,line_4=line_4,img_4=img_4,conv=conv, wayemi=wayemi,imgfront=imgfront,filenm=filenm,gotdir=gotdir,file_abund=file_abund,vers=vers'
   return
endif

; INPUT:
; rho: (2d float array) number density in CGS
; tem: (2d float array) temperature in CGS 
; v1m, v2m: (2d float arrays) x and y components of velocity, in km/s 
; gridx, gridy: (floats) x and y axes 
; ion: (string) acronym of the ion
; w0: (float) wavelength of line center 
; mua_d: (floats) array of line-of-sight angles. The routine is set to return
;        a maximum of 4 angles. If more are necessary you just need to add
;        the corresponding keywords in the call (just changing the number).
; gotdir: (string) directory path to where the .dat G(T,n) file is. 
; file_abund: (string) file for abundance abundance. 2 kinds are implemented:
;            'photospheric' or 'coronal' corresponding, respectively, to the
;            CHIANTI packages: sun_coronal.abund and sun_photospheric.abund
; vers: (int) the CHIANTI version (6 or 7).
; wayemi: (int) for different ways of calculating the emissivity in the
;        routine lineongrid_goft_tab.pro. The default is wayemi =
;        4. Set wayemi = 5 when only imaging is desired (for SDO/AIA
;        filters for instance)

; Optional:
; conv: set for converting density into number density when the latter is
;       in SI units, has been normalized by 1.e10 and the plasma is
;       fully ionized
; filenm = (string) relevant only for imaging purposes (wayemi = 5). Used for
;         saving the emission_goft file in routine interpol_emiss_data.pro

; OUTPUT:
; gridz & imgfront: (1d and 2D float arrays, resp.) grid along the z-axis. This is only relevant when the
;        parameter 'width' defined below is greater than 0. By doing so you
;        define a depth for your 2D plane, and the parameter 'imgfront' will
;        therefore represent an intensity image of the 2D plane. 
; emission_goft: (2d float array) calculated emissivity values at each point point of the 2D
;        plane. 
; wave: (1d float array) the wavelength array for the line transition of interest
; nwave: (float) number of points in wave array as set in routine
;       lineongrid_goft_tab.pro. nwave is 100 points by default
; n_gridx_<num>: (1d float array) x-axis coordinates for the new grid corresponding to
;           line-of-sight angle <num>. This is an array with a specific order
;           produced by the routine gridlos.pro.
; n_gridy_<num>: same as above for the y-coordinates.
; ngrid_<num>: (1d int array) the i-th position of this array contains the number of
;             points of the i-th ray for line_of_sight <num>. See
;             gridlos.pro for more details.
; dl_<num>: (float) the spatial resolution of the grid for line-of-sight
;          <num>, as set by gridlos.pro
; line_1: 2D float array containing the line profile for each line-of-sight
;        ray, obtained after integrating along the line-of-sight (calls
;        routine integrateemission.pro). The dimensions are thus (number of
;        rays, number of wavelength points) = (n_elements(ngrid_<num>), nwave)

; CALLS:
; lineongrid_goft_tab, integrateemission, interpol_emiss_data, gridlos

if keyword_set(conv) then begin
   proton = 1.67262158*10^(-27.)
   kboltz = 1.380658*10^(-23.)
   normro = 1.e10
   normte = proton/(2*kboltz)*normro
;   nefrac = 1.2
   nefrac = 1.
   rh = rho*normro*nefrac
   te = tem /normte
;  ne_s = rh/proton/1.e6 
   nem = ne_s
endif else begin
   ; check for CGS
   rh = rho
   te = tem
   ne_s = nem
endelse

direction = 4
imaging = 0
nang = n_elements(mua_d)
width = 0.

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
      if mua eq 0. then direction = 1
      if mua eq 90. then direction = 2
      if mua ne 0. and mua ne 90. then direction = 4
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

      void = execute(exe1)
      void = execute(exe2)
      void = execute(exe3)
      void = execute(exe4)
      void = execute(exe5)
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
      if mua eq 0. then direction = 1
      if mua eq 90. then direction = 2
      if mua ne 0. and mua ne 90. then direction = 4
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

      void = execute(exe1)
      void = execute(exe2)
      void = execute(exe3)
      void = execute(exe4)
      void = execute(exe5)
   endfor
endelse

end
