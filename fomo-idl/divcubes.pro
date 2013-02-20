

pro divcubes,dimt=dimt,dimx=dimx,dimy=dimy,model=model,sngcub=sngcub,mag=mag,smooth=smooth

if keyword_set(model) eq 0 then begin
   print,'Check input and output directories first'
   print,'divcubes,dimt=dimt,dimx=dimx,dimy=dimy,model=model[,mag=mag,smooth=smooth]'
   return
endif

; Divides cubes produced by datacubes_wt.pro into slices of cylinder
; for each x position (up to half the cylinder) and time step
; INPUT:
; dimt = number of points in time dimension
; dimx = dimension in x direction
; dimy = dimension in y direction
; model = string with name of treated model:
;      model = 'base'corresponds to ka = 2.24, 
;      model = 'long' corresponds to ka = 1.25
;      model = 'high_res' corresponds to ka = 2.24, high spatial resolution
;      model = 'high_res2d' corresponds to ka = 2.24, high 2D spatial resolution
;      model = 'highT' corresponds to ka = 2.24, high external temperature
; OUPUT:
; Produces IDL save files of slices at each x and time position
; OPTIONAL:
; if keyword 'mag' is set then it includes the magnetic field cubes
; if keyword 'smooth' is set then cubes correspond to smooth profiles

  if keyword_set(smooth) then sm = 1 else sm = 0
  if model eq 'base' then  kanm = 'ka2.24'
  if model eq 'long' then  kanm = 'ka1.25'
  if model eq 'high_res' then  kanm = 'ka2.24_hgres'
  if model eq 'high_res2d' then  kanm = 'ka2.24_hgres2d'
  if model eq 'highT' then kanm = 'ka2.24highT'

; INPUT DIRECTORY:
  cubedir='../cubes/set2/'
; OUTPUT DIRECTORY:
  savedir='../cubes/set2/'
;  dir = '/volume1/scratch/set3/'

  for i=0,dimt-1 do begin
     restore,cubedir+'cubes_'+kanm+'_'+string(i,format='(i3.3)')+'.sav'
     if sm eq 0 then begin
        for j=0,dimx/2 do begin
           if sngcub eq 'all' or sngcub eq 'te' then te = te_cube[j,*,*]
           if sngcub eq 'all' or sngcub eq 'rh' then rho = rh_cube[j,*,*]
           if sngcub eq 'all' or sngcub eq 'vr' then vr = vr_cube[j,*,*]
           if sngcub eq 'all' or sngcub eq 'vz' then vz = vz_cube[j,*,*]
           if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bz')) then bz = bz_cube[j,*,*]
           if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'br')) then br = br_cube[j,*,*]
           save,te,rho,vr,vz,br,bz,filename=savedir+'slice_'+sngcub+'_'+kanm+'_'+string(i,format='(i3.3)')+'t_'+string(j,format='(i3.3)')+'x'+'.sav'
        endfor
     endif
     if sm eq 1 then begin
        for j=0,dimx/2 do begin
           if sngcub eq 'all' or sngcub eq 'te' then te = te_cube_sm[j,*,*]
           if sngcub eq 'all' or sngcub eq 'rh' then rho = rh_cube_sm[j,*,*]
           if sngcub eq 'all' or sngcub eq 'vr' then vr = vr_cube_sm[j,*,*]
           if sngcub eq 'all' or sngcub eq 'vz' then vz = vz_cube_sm[j,*,*]
           if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'bz')) then bz = bz_cube_sm[j,*,*]
           if (keyword_set(mag) and (sngcub eq 'all' or sngcub eq 'br')) then br = br_cube_sm[j,*,*]
           save,te,rho,vr,vz,br,bz,filename=savedir+'slice_'+sngcub+'_'+kanm+'_'+string(i,format='(i3.3)')+'t_'+string(j,format='(i3.3)')+'x'+'.sav'
        endfor
     endif
     print,string(13b)+' % finished: ',float(i)*100./(dimt-1),format='(a,f4.0,$)'
  endfor

end
