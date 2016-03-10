
pro wavetube_synth,dir=dir,gotdir=gotdir,w0=w0,ion=ion,idcase=idcase,sl=sl,imaging=imaging,aia=aia,save=save

; WRAPPER AROUND FOMO:
; This example assumes that you have 2 IDL save files (wavetube0.idl &
; wavetube1.idl) with your numerical data (CGS UNITS), correpsonding to two
; different time steps. The data is organized as [x,y,z] and the 
; forward modelling is done by taking slices in the (x,y) plane. 
; The angles provided below correspond to angles in this plane, 
; where 0 degrees corresponds to the x-axis, 90 degrees to the 
; y-axis (positive).

; Numerical data needed: number density, temperature and velocity
; components in the LOS planes.

; KEYWORDS:
; DIR: work directory. (don't forget to add '/' at the end of the path
; IDCASE: name of the current study. The wrapper assumes that a
; subdirectory in dir exists with this name.
; SL: subdirectory indicating the slice orientation (ex: sl='2dxy')
; The output will be saved in dir+idcase+'/'+sl+'/sav/'
; W0 = line center for emission line of interest
; ION = name of line (see below for a partial list, full list in FoMo/fomo-idl/elements.pro)
; IMAGING: set this keyword indicating that only imaging (no spectra)
; will be saved
; AIA: keyword indicating that AIA channels will be treated. See below
; for details on how to select the appropriate channel
; SAVE: Set this keyword to save all output

if keyword_set(dir) eq 0 then begin
    print,'wavetube_synth,dir=dir,gotdir=gotdir,w0=w0,ion=ion,idcase=idcase,sl=sl,imaging=imaging,aia=aia,save=save'
    return
endif

; NAME of output:
if keyword_set(imaging) then name = 'wavetube_imag_'+idcase+'_'+sl else name = 'wavetube_synth_'+idcase+'_'+sl

; OUTPUT DIRECTORY:

savedir = dir+idcase+'/'+sl+'/sav/'
if (not file_test(savedir,/directory)) then begin
   print,'savedir does not exist. Please check path.'
   return
endif

; CHIANTITABLES:
if (not file_test(gotdir,/directory)) then begin
   print,'Chiantitables directory does not exist. Please check path.'
   return
endif

; ABUNDANCE: 'coronal' or 'photospheric'
file_abund = 'coronal'

; CHOOSE LINE TO SYNTHESIZE:

; for Fe IX 171.073
;w0 = 171.073 & ion = 'fe_9'

; for Fe XII 193.509
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
;w0 = 193. & ion = '193'
; AIA-193
;w0 = 193. & ion = '193'
; AIA-211
;w0 = 211. & ion = '211'
; AIA-131
;w0 = 131. & ion = '131'
; AIA-304
;w0 = 304. & ion = '304'
; AIA-335
;w0 = 335. & ion = '335'
; AIA-094
;w0 = 94. & ion = '094'

nw0 = string(w0,format='(d0.3)')

; CHOOSE ANGLES:
; 0-deg corresponds to x-axis, 90-deg to y-axis
; up to 7 angles accepted
; mua_d = [0.,15.,30.,45.,60.,75.,90.]

mua_d = [0.,45.,90.]

; SDO FILTERS (imaging): 
if keyword_set(aia) or keyword_set(imaging) then begin
   if keyword_set(aia) then begin
      bilinear = 1
      wayemi = 5
      name = name+'aia_'+ion
   endif
   if keyword_set(imaging) then begin
      wayemi = 5
      name = name +'_'+ion+'_'+nw0
   endif
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
inin3 = 0 ; initial slice
n3 = 256 ; total number

; Go through all, indicating which one the current idl thread is
; working on:

for i=0,num-1 do begin
   if num lt 1000 then inum = STRING(i, FORMAT = "(I3.3)") else inum = STRING(i, FORMAT = "(I4.4)")
   sn = string(i,format="(i1)")

; location of your save files with numerical data. In this case it
; reads a 3D cube [512x512x256] corresponding to [x,y,z] for each time
; step:
   restore,dir+idcase+'/wavetube'+sn+'.idl'
   gridx = x ; x-grid
   gridy = y ; y-grid
   for j=inin3,n3-1 do begin
      jnum = STRING(j, FORMAT = "(I3.3)")
      lcx3 = j
      if n3 gt 1 then checkfile = name+'_'+inum+'_'+jnum+'.done' else checkfile = name+'_'+inum+'.done'
      if (file_test(savedir+checkfile) eq 0) then begin
         OPENW, lun, savedir+checkfile, /get_lun
         free_lun, lun
         print,'Doing slice_time:',checkfile
         nem = ndens[*,*,j]
         tem = temperature[*,*,j]
         v1m = v_x[*,*,j]/1.e5
         v2m = v_y[*,*,j]/1.e5

         ; CALCULATE INTENSITIES
         synthemi,rho=rho,nem=nem,tem=tem,v1m=v1m,v2m=v2m,ion=ion,mua_d=mua_d,gridx=gridx,gridy=gridy,emission_goft=emission_goft,wave=wave,nwave=nwave,w0=w0,n_gridx_1=n_gridx_1,n_gridy_1=n_gridy_1,ngrid_1=ngrid_1,n_gridx_2=n_gridx_2,n_gridy_2=n_gridy_2,ngrid_2=ngrid_2,n_gridx_3=n_gridx_3,n_gridy_3=n_gridy_3,ngrid_3=ngrid_3,n_gridx_4=n_gridx_4,n_gridy_4=n_gridy_4,ngrid_4=ngrid_4,n_gridx_5=n_gridx_5,n_gridy_5=n_gridy_5,ngrid_5=ngrid_5,n_gridx_6=n_gridx_6,n_gridy_6=n_gridy_6,ngrid_6=ngrid_6,n_gridx_7=n_gridx_7,n_gridy_7=n_gridy_7,ngrid_7=ngrid_7,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,dl_4=dl_4,dl_5=dl_5,dl_6=dl_6,dl_7=dl_7,line_1=line_1,line_2=line_2,line_3=line_3,line_4=line_4,line_5=line_5,line_6=line_6,line_7=line_7,wayemi=wayemi,filenm=filenm,gotdir=gotdir,file_abund=file_abund

         if keyword_set(save) then begin
            if i eq 0 and j eq inin3 then save,gridx,gridy,wave,w0,ion,mua_d,filename=savedir+'params_'+name+'.sav'

            save,time,dl_1,dl_2,dl_3,dl_4,dl_5,dl_6,dl_7,ngrid_1,ngrid_2,ngrid_3,ngrid_4,ngrid_5,ngrid_6,ngrid_7,n_gridx_1,n_gridx_2,n_gridx_3,n_gridx_4,n_gridx_5,n_gridx_6,n_gridx_7,n_gridy_1,n_gridy_2,n_gridy_3,n_gridy_4,n_gridy_5,n_gridy_6,n_gridy_7,line_1,line_2,line_3,line_4,line_5,line_6,line_7,emission_goft,filename=savedir+name+'_'+inum+'_'+jnum+'.sav'
         endif
         print,string(13b)+' % finished: ',float((i+1)*j)*100./(num*n3-1),format='(a,f4.0,$)'
      endif
   endfor

endfor

end

