
pro fomo_synth,dir=dir,gotdir=gotdir,w0=w0,ion=ion,idcase=idcase,imaging=imaging,channel=channel,filenm=filenm,extro=extro,silent=silent,save=save

if keyword_set(dir) eq 0 then begin
   print,'fomo_synth,dir=dir,gotdir=gotdir,w0=w0,ion=ion,idcase=idcase,imaging=imaging,channel=channel,filenm=filenm,extro=extro,silent=silent,save=save'
   return
endif

; WRAPPER AROUND FOMO:
; This example assumes that you have 2 IDL save files (wavetube0.idl &
; wavetube1.idl) with your numerical data (CGS UNITS), correpsonding to two
; different time steps. The data is organized as [x,y,z] and the 
; forward modelling is done by taking slices in the (x,y)
; plane. Specifically wavetube contains (all in CGS and all have
; dimensions (dimx,dimy,dimz) ):
; - ndens: total number density (converted to electron number density
;   during calculation, assuming ratio of protons to electrons is
;   0.887 and 0.848, respectively, for coronal and photospheric abundances)
; - temperature
; 
; In this example we further assume that we have 50 positions along z
; (leading to 50 slices for each time step)
; The angles provided below correspond to angles in this plane, 
; where 0 degrees corresponds to the x-axis, 90 degrees to the 
; y-axis (positive).

; Numerical data needed: number density, temperature and velocity
; components in the LOS planes.

; KEYWORDS:
; DIR: work directory. (don't forget to add '/' at the end of the path
; GOTDIR: directory of chianti tables in fomo (full path)
; IDCASE: name of the current study. The wrapper assumes that a
; subdirectory in dir exists with this name.
; The output will be saved in dir+idcase+'/sav/'
; W0 = line center for emission line of interest
; ION = name of line (see below for examples) 
; IMAGING: set this keyword indicating that only imaging (no spectra)
; will be saved
; CHANNEL: set to 'aia' for AIA channels. See below
; for details on how to select the appropriate channel
; FILENM = optional keyword for additional labelling
; EXTRO: set for G(T,n) tables with extended density range [6,12]
;        in log. Default for spectral lines is [8,11] for log(T)>5 and [8,12]
;        for log(T)<5 where T is the maximum formation temperature. 
;        Default for AIA 304,1600,1700,4500 is [8,12] and [8,11] for rest.
;        Default for EIT 304 is [8,12] and [8,11] for rest.
;        Default for DKIST is [8,11].
; SAVE: Set this keyword to save all output

if keyword_set(dir) eq 0 then begin
    print,'fomo_synth,dir=dir,gotdir=gotdir,w0=w0,ion=ion,idcase=idcase,imaging=imaging,channel=channel,filenm=filenm,silent=silent,save=save'
    return
endif
if keyword_set(channel) then imaging = 1
; NAME of output:
if (keyword_set(imaging) or keyword_set(channel)) then name = 'fomo_'+idcase+'_imag' else name = 'fomo_'+idcase+'_synth'

; OUTPUT DIRECTORY:

savedir = dir+idcase+'/sav/'
file_mkdir,savedir
print,'Save directory: '+savedir

; CHIANTITABLES:
if (not file_test(gotdir,/directory)) then begin
   print,'Chiantitables directory does not exist. Please check path.'
   return
endif

; ABUNDANCE: 'coronal' or 'photospheric'
file_abund = 'coronal'

; CHOOSE LINE TO SYNTHESIZE:
; examples:
; for Fe IX 171.073 (about 1MK)
;w0 = 171.073 & ion = 'fe_9'

; for Fe XII 193.509 (about 1.5MK)
;w0 = 193.509 & ion = 'fe_12'

; Si VII 197.7684
;w0 = 197.768 & ion = 'si_7'

; He II 303.781
;w0 = 303.781 & ion = 'he_2'

; He II 303.786
;w0 = 303.786 & ion = 'he_2'

; C IV 1548.189
;w0 = 1548.1899 & ion = 'c_4'

; C IV 1550.775
;w0 = 1550.775 & ion = 'c_4'

; Ne VIII 770.4103
;w0 = 770.4103 & ion = 'ne_8'

; O IV  1399.78
;w0 = 1399.78 & ion = 'o_4'

; O IV  1401.16
;w0 = 1401.16 & ion = 'o_4'

; Si IV 1393.757
;w0 = 1393.757 & ion = 'si_4'

; Si IV 1402.772
;w0 = 1402.772 & ion = 'si_4'

; C II 1334.535
;w0 = 1334.535 & ion = 'c_2'

; C II 1335.71
;w0 = 1335.71 & ion = 'c_2'

; Mg II h & k: 2796.3521
;w0 = 2796.3521 & ion = 'mg_2'

; Mg II h & k: 2803.531
;w0 = 2803.531 & ion = 'mg_2'


; SDO FILTERS:
; AIA-171
; w0 = 171
; AIA-193
; w0 = 193
; AIA-211
; w0 = 211
; AIA-131
; w0 = 131
; AIA-304
; w0 = 304
; AIA-335
; w0 = 335
; AIA-094
; w0 = 94
; AIA-1600
; w0 = 1600
; AIA-1700
; w0 = 1700

if ~keyword_set(ion) and keyword_set(channel) and keyword_set(w0) then begin
   if w0 lt 1.e3 then ion = string(w0,format="(i3.3)") else ion = string(w0,format="(i4.4)")
endif

nw0 = string(w0,format='(d0.3)')

; CHOOSE ANGLES:
; 0-deg corresponds to x-axis, 90-deg to y-axis
; up to 7 angles accepted
; mua_d = [0.,15.,30.,45.,60.,75.,90.]

mua_d = [0.,45.,90.]
nang = n_elements(mua_d)
for i=0,nang-1 do begin
   sn = string(i+1,format="(i1)")
   if i eq 0 then linlin = 'line_'+sn else linlin = linlin+',line_'+sn
   if i eq 0 then dllin = 'dl_'+sn else dllin = dllin+',dl_'+sn
   if i eq 0 then nglin = 'ngrid_'+sn else nglin = nglin+',ngrid_'+sn
   if i eq 0 then ngxlin = 'n_gridx_'+sn else ngxlin = ngxlin+',n_gridx_'+sn
   if i eq 0 then ngylin = 'n_gridy_'+sn else ngylin = ngylin+',n_gridy_'+sn
endfor
emilin = 'emission_goft'
savlin = linlin+','+emilin

; SDO FILTERS (imaging): 
if keyword_set(channel) or keyword_set(imaging) then begin
   wayemi = 5
   if keyword_set(channel) then begin
      name = name+'_'+channel+'_'+ion
   endif else begin
      name = name +'_'+ion+'_'+nw0
   endelse
endif else begin
   wayemi = 4
   name = name +'_'+ion+'_'+nw0
endelse

;UNITS:
; make sure you have CGS units for all

;proton=1.67262158*10^(-24.)
;kboltz = 1.380658*10^(-16.)
;gamma = 5./3.

; TIMESTEPS and NUMBER OF SLICES:

; Indicate the number of time steps:
num = 2

; GRID:
; make sure you have uniform grids

; Indicate the number of slices (here, along z-axis)
nz0 = 0 ; initial slice
dimz = 50 ; total number

; Go through all, indicating which one the current idl thread is
; working on:

for i=0,num-1 do begin
   sn = string(i,format="(i1)")
   inm = STRING(i, FORMAT = "(I4.4)")
; location of your save files with numerical data. In this case it
; reads a 3D cube [512x512x50] corresponding to [x,y,z] for each time
; step:

   restore,dir+idcase+'/wavetube'+sn+'.idl'
   gridx = x ; x-grid
   gridy = y ; y-grid
   for j=nz0,dimz-1 do begin
      jnm = STRING(j, FORMAT = "(I5.5)")
      lcx3 = j
      checkfile = name+'_z='+jnm+'_t='+inm+'.done'
      if (file_test(savedir+checkfile) eq 0) then begin
         OPENW, lun, savedir+checkfile, /get_lun
         free_lun, lun
         print,'Doing slice_time:',checkfile
; ndens, temperature, v_x, v_y are 3d-arrays (dimx,dimy,dimz)
         nem = ndens[*,*,j] ; total density (CGS)
         tem = temperature[*,*,j] ; in K
         if ~keyword_set(imaging) and ~keyword_set(channel) then begin
            v1m = v_x[*,*,j]/1.e5 ; convert to km/s
            v2m = v_y[*,*,j]/1.e5 ; convert to km/s
         endif
         tstep = t ; t is scalar indicating time step in sec

         ; CALCULATE INTENSITIES
         synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,n_gridx_5=n_gridx_5,n_gridy_5=n_gridy_5,ngrid_5=ngrid_5,n_gridx_6=n_gridx_6,n_gridy_6=n_gridy_6,ngrid_6=ngrid_6,n_gridx_7=n_gridx_7,n_gridy_7=n_gridy_7,ngrid_7=ngrid_7,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,dl_5=dl_5,dl_6=dl_6,dl_7=dl_7,line_1=line_1,line_2=line_2,line_3=line_3,line_4=line_4,line_5=line_5,line_6=line_6,line_7=line_7,wayemi=wayemi,filenm=filenm,gotdir=gotdir,file_abund=file_abund,imaging=imaging,channel=channel,nab=nab,abund_name=abund_name,enum=enum,inum=inum,abund_fact=abund_fact,extro=extro,silent=silent

         if keyword_set(save) then begin
            if i eq 0 and j eq nz0 then begin
               exp = 'save,gridx,gridy,wave,w0,ion,mua_d,dimz,num,'+dllin+','+nglin+','+ngxlin+','+ngylin+',filename=savedir+"params_"+name+".sav"'
               void = execute(exp)
            endif
            exs = 'save,tstep,'+savlin+',filename=savedir+name+"_z="+jnm+"_t="+inm+".sav"'
            void = execute(exs)
         endif
      endif
   endfor
   print,string(13b)+' % finished: ',i*100./(float(num)-1),format='(a,f4.0,$)'
endfor

end

