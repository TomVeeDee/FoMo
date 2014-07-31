pro buildcubes

; Builds data cubes from interpolated slices in C++

dimx=180.
dimy=180.
dimz=210.
nvar = 5.
harmonic=3
dimt= 48./harmonic

norm = 1.e6
mup = 1.25663706*1.e-6/norm
kboltz = 1.380658*10^(-23.)
proton = 1.67262158*10^(-27.)

;t=1
for t=0,dimt-1 do begin
  
  rho = fltarr(dimx,dimy,dimz)
  te = fltarr(dimx,dimy,dimz)
  br = fltarr(dimx,dimy,dimz)
  bt = fltarr(dimx,dimy,dimz)
  bz = fltarr(dimx,dimy,dimz)

  for z=0, dimz-1 do begin
    
    dataslice = fltarr(2+nvar,dimx,dimy)

    eigfcube = '/users/cpa/sgijsen/idl_general/test_advectioncube/eigftcube'+string(t,format="(i3.3)")+'z'+string(z,format="(i3.3)")+'.dat'
    openr,lun,eigfcube,/get_lun
    readf,lun,dataslice
    close,lun
    free_lun,lun
    
    rho[*,*,z] = transpose(reform(dataslice[2,*,*]))/(1.e6 * proton *norm^3)
    te[*,*,z] = transpose(reform(dataslice[3,*,*]))*proton/kboltz*norm^2*3./8.
    br[*,*,z] = transpose(reform(dataslice[4,*,*]))
    bt[*,*,z] = transpose(reform(dataslice[5,*,*]))
    bz[*,*,z] = transpose(reform(dataslice[6,*,*]))
    
  endfor
  print, 'time step '+strtrim(string(t),1)
  save, rho, te, br, bt, bz, filename='/users/cpa/sgijsen/adv/razin/eigfhar'+strtrim(string(harmonic),1)+'t'+string(t,format="(i3.3)")+'.sav'
  
endfor
end