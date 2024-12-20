
pro synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,gridz=gridz,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,n_gridx_5=n_gridx_5,n_gridy_5=n_gridy_5,ngrid_5=ngrid_5,n_gridx_6=n_gridx_6,n_gridy_6=n_gridy_6,ngrid_6=ngrid_6,n_gridx_7=n_gridx_7,n_gridy_7=n_gridy_7,ngrid_7=ngrid_7,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,dl_5=dl_5,dl_6=dl_6,dl_7=dl_7,line_1=line_1,line_2=line_2,line_3=line_3,line_4=line_4,line_5=line_5,line_6=line_6,line_7=line_7,wayemi=wayemi,filenm=filenm,gotdir=gotdir,file_abund=file_abund,ext_abund=ext_abund,imaging=imaging,revvel=revvel,channel=channel,nab=nab,abund_name=abund_name,enum=enum,inum=inum,abund_fact=abund_fact,extro=extro,silent=silent

if arg_present(rho) lt 1 and arg_present(nem) lt 1 then begin
   print,'synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,n_gridx_5=n_gridx_5,n_gridy_5=n_gridy_5,ngrid_5=ngrid_5,n_gridx_6=n_gridx_6,n_gridy_6=n_gridy_6,ngrid_6=ngrid_6,n_gridx_7=n_gridx_7,n_gridy_7=n_gridy_7,ngrid_7=ngrid_7,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,dl_5=dl_5,dl_6=dl_6,dl_7=dl_7,line_1=line_1,line_2=line_2,line_3=line_3,line_4=line_4,line_5=line_5,line_6=line_6,line_7=line_7,wayemi=wayemi,filenm=filenm,gotdir=gotdir,file_abund=file_abund,ext_abund=ext_abund,imaging=imaging,revvel=revvel,channel=channel,nab=nab,abund_name=abund_name,enum=enum,inum=inum,abund_fact=abund_fact,extro=extro,silent=silent'
   return
endif


; INPUT:
; rho (or nem): (2d float array) total density (or total number density) in CGS
; tem: (2d float array) temperature in CGS 
; v1m, v2m: (2d float arrays) x and y components of velocity, in km/s 
; gridx, gridy: (floats) x and y axes in cgs, corresponding to plane
; of slice passed to fomo. Note: they are assumed uniform
; gridz: (float) z axis in cgs, corresponding to perpendicular plane
; to slice passed to fomo (assumed uniform).
; ion: (string) acronym of the ion
; w0: (float) wavelength of line center 
; mua_d: (floats) array of line-of-sight angles. The routine is set to return
;        a maximum of 4 angles. If more are necessary you just need to add
;        the corresponding keywords in the call (just changing the number).
; gotdir: (string) directory path to where the .dat G(T,n) file is. 
; file_abund: (string) file for abundance abundance. 2 kinds are implemented:
;            'photospheric' or 'coronal' corresponding, respectively, to the
;            CHIANTI packages: sun_coronal.abund and sun_photospheric.abund
;            By default the 'coronal' abundance package is set.
;	     If other abundance is desired set file_abund = 'other' and provide 
;	     the full path to the abundance file in the keyword 'ext_abund'.
; imaging: set this keyword for only imaging (no spectral info)
; wayemi: (int) for different ways of calculating the emissivity in the
;        routine lineongrid_goft_tab.pro. The default is wayemi =
;        4. Set wayemi = 5 when only imaging is desired (for SDO/AIA
;        filters for instance)
; revvel: By default the LOS angles are set considering the
; x-direction and y-directions as positive. If this is not the case
; then the velocity is of opposite signs with respect to the LOS. Set
; this keyword in this case.
; channel: imaging channel. e.g. channel='aia' for SDO/AIA filters (imaging is then set to 1 automatically)
; enum = nuclear charge of element
; inum = ionisation stage of element
; nab = indicates abundance package '_abph', '_abco' or '_abext'
; abund_name = full path to selected abundance package
; abund_fact = factor for multiplication with abundance

; Optional:

; filenm = (string) relevant only for imaging purposes (wayemi = 5). Used for
;         saving the emission_goft file in routine interpol_emiss_data.pro
;         Also used for specific labelling of chianti tables.
; extro: set for G(T,n) tables with extended density range [6,12]
;        in log. Default for spectral lines is [8,11] for log(T)>5 and [8,12]
;        for log(T)<5 where T is the maximum formation temperature. 
;        Default for AIA 304,1600,1700,4500 is [8,12] and [8,11] for rest.
;        Default for EIT 304 is [8,12] and [8,11] for rest.
;        Default for DKIST is [8,11].

; OUTPUT:
 
; emission_goft: (2d float array) calculated emissivity values at each
;        point point of the 2D plane. 
; wave: (1d float array) the wavelength array for the line transition of interest
; nwave: (float) number of points in wave array as set in routine
;       lineongrid_goft_tab.pro. nwave is 100 points by default
; n_gridx_<num>: (1d float array) x-axis coordinates for the new grid corresponding to
;           line-of-sight angle <num>. This is an array with a specific order
;           produced by the routine gridlos.pro.
; n_gridy_<num>: same as above for the y-coordinates.
; ngrid_<num>: (1d int array) the i-th position of this array contains
;        the number of points of the i-th ray for line_of_sight
;        <num>. See gridlos.pro for more details.
; dl_<num>: (float) the spatial resolution of the grid for line-of-sight
;          <num>, as set by gridlos.pro
; line_1: 2D float array containing the line profile for each line-of-sight
;        ray, obtained after integrating along the line-of-sight (calls
;        routine integrateemission.pro). The dimensions are thus (number of
;        rays, number of wavelength points) = (n_elements(ngrid_<num>), nwave)

; CALLS:
; lineongrid_goft_tab, integrateemission, interpol_emiss_data, gridlos

proton = 1.67262158*10^(-24.)

if keyword_set(file_abund) then begin
   if file_abund eq 'photospheric' then fnhne = 0.848 ; ratio of protons to electrons number (= proton_dens(6.0) from Chianti, with photospheric abundances)
   if file_abund eq 'coronal' then fnhne = 0.887 ; ratio of protons to electrons number (= proton_dens(6.0) from Chianti, with coronal abundances)
   if file_abund eq 'other' then begin
      if ~keyword_set(ext_abund) then begin
         print,'Please provide the name and full path of the abundance file in the keyword "ext_abund". Note: ratio of protons to electrons taken as 0.887. Modify accordingly.'
         stop
      endif else begin
         abund_name = ext_abund
         nab = '_abext'
         fnhne = 0.887
      endelse
   endif
endif else begin
   fnhne = 0.887 ; ratio of protons to electrons number (= proton_dens(6.0) from Chianti, with coronal abundances)
endelse

te_s = tem
if ~keyword_set(nem) and keyword_set(rho) then ne_s = rho/(proton*fnhne)
if keyword_set(nem) then ne_s = nem/(1.+fnhne)

direction = 4
nang = n_elements(mua_d)

if ~keyword_set(imaging) and ~keyword_set(channel) then begin
   velx = v1m
   vely = v2m
endif
ngridsx = n_elements(gridx)
ngridsy = n_elements(gridy)
dx = gridx[1]-gridx[0]
dy = gridy[1]-gridy[0]
dz = gridz[1]-gridz[0]
d_perp = dz

if keyword_set(imaging) or keyword_set(channel) then begin
   ; Channel imaging (filters) or line imaging:
   wayemi = 5
   logt = alog10(tem)
   if keyword_set(channel) then interpol_emiss_data,ne_s=ne_s,te=tem,ion=ion, w0=w0,emission_goft=emission_goft,filenm=filenm,file_abund=file_abund,gotdir=gotdir,channel=channel,extro=extro,silent=silent else lineongrid_goft_tab, te_s=te_s, ne_s=ne_s, gotdir=gotdir,wave=wave,nwave=nwave,ion=ion, w0=w0, emission_goft=emission_goft,goft=goft, logt=logt,wayemi=wayemi,watom=watom,file_abund=file_abund,ext_abund=ext_abund,nab=nab,abund_name=abund_name,enum=enum,inum=inum,abund_fact=abund_fact,extro=extro,silent=silent
   
   for i=0,nang-1 do begin
      mua = mua_d[i]
      if mua eq 0. then direction = 2
      if mua eq 90. then direction = 1
      if mua ne 0. and mua ne 90. then direction = 4
      gridlos, gridx=gridx, gridy=gridy, mua_d=mua, dx=dx, dy=dy, n_gridx=n_gridx, n_gridy=n_gridy, ngrid=ngrid, dl=dl,ds=ds

      integrateemission,emission=emission_goft,logt=logt,n_gridx=n_gridx,n_gridy=n_gridy,ngrid=ngrid,w0=w0,direction=direction,imaging=imaging,channel=channel,image=image,dl=dl,ds=ds,d_perp=d_perp,wayemi=wayemi,silent=silent
      dlos = (size(image))[1]
      inan = string(i+1,format="(i1)")
      exe1 = 'line_'+inan+' = image'
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

   lineongrid_goft_tab, te_s=te_s, ne_s=ne_s, gotdir=gotdir,wave=wave,nwave=nwave,ion=ion, w0=w0, emission_goft=emission_goft,goft=goft, logt=logt,wayemi=wayemi,watom=watom,file_abund=file_abund,ext_abund=ext_abund,nab=nab,abund_name=abund_name,enum=enum,inum=inum,abund_fact=abund_fact,extro=extro,silent=silent

   for i=0,nang-1 do begin
      mua = mua_d[i]
      if mua eq 0. then direction = 2
      if mua eq 90. then direction = 1
      if mua ne 0. and mua ne 90. then direction = 4
      gridlos, gridx=gridx, gridy=gridy, mua_d=mua, velx=velx, vely=vely, dx=dx, dy=dy, n_gridx=n_gridx, n_gridy=n_gridy, ngrid=ngrid, dl=dl,ds=ds, losvel=losvel

      if keyword_set(revvel) then begin
         losvel = -losvel
      endif else begin 
         losvel = losvel
      endelse

      integrateemission,emission=emission_goft,logt=logt,n_gridx=n_gridx,n_gridy=n_gridy,ngrid=ngrid,wave=wave,w0=w0,direction=direction,losvel=losvel,imaging=imaging,image=image,watom=watom,wayemi=wayemi,dl=dl,ds=ds,d_perp=d_perp,silent=silent

      dlos = (size(image))[1]
      inan = string(i+1,format="(i1)")
      exe1 = 'line_'+inan+' = image'
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
