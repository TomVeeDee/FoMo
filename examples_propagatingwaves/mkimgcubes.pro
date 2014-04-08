
pro mkimgcubes, set=set, gridx=gridx, dimt=dimt, dim00d=dim00d, dim30d=dim30d, dim45d=dim45d, dim60d=dim60d, gridy, gridz_ext, wave, imgcube00d=imgcube00d,imgcube30d=imgcube30d,imgcube45d=imgcube45d,imgcube60d=imgcube60d,imaging=imaging

if keyword_set(set) eq 0 then begin
   print,'mkimgcubes, set=set, gridx=gridx, dimt=dimt, dim00d=dim00d, dim30d=dim30d, dim45d=dim45d, dim60d=dim60d, gridy, gridz_ext, wave, imgcube00d,imgcube30d,imgcube45d,imgcube60d'
   return
endif

; Produces 3D intensity image cubes from sav files produced by prl_slices.pro
; INPUT:
; set = defines model:
;     = 2: 'base model - 171' corresponds to ka = 2.24, line Fe IX 171
;     = 21: 'high T model' base model with high external temperature
;     = 22: 'smooth model' base model with smooth density profile across cylinder
;     = 23: 'base model - 193' base model with line Fe XII 193
;     = 3 : 'long lambda model' corresponds to ka = 1.25
;     = 24: 'high_res' corresponds to base model with high resolution
;     = 25: 'high_res2d.fe9': 2D base model with high resolution, 171 emission
;     = 26: 'high_res2d.fe12' 2D base model with high resolution, 193 emission
; gridx = grid along x axis (radial direction)
; dimt = number of points in time dimension
; dim00d = number of points in longitudinal direction in new grid with angle=00
; dim30d = number of points in longitudinal direction in new grid with angle=30
; dim45d = number of points in longitudinal direction in new grid with angle=45
; dim60d = number of points in longitudinal direction in new grid with angle=60
; set keyword imaging if no spectral info

if keyword_set(imaging) then imaging = 1 else imaging = 0

  savcub = 1
  if set ne 25 and set ne 26 then dimx = n_elements(gridx) else dimx = 1
;  if set eq 2 then begin
;     dim00d = 342 & dim30d = 398 & dim45d = 386 & dim60d = 348
;  endif
;  if set eq 3 or set eq 31 then begin
;     dim00d = 404 & dim30d = 452 & dim45d = 430 & dim60d = 379
;  endif


if set eq 2 then begin dir = '/volume1/scratch/set2/' & ka_n = '_ka2.24_fe9' & sln = 3 & endif
;if set eq 2 then begin dir = '/users/cpa/pantolin/Modeling/cubes/set2/' & ka_n = '_ka2.24_fe9' & sln = 6 & endif
;if set eq 21 then begin dir = '/users/cpa/pantolin/Modeling/cubes/set2/highT/' & ka_n = '_ka2.24_highT_fe9' & sln = 7 & endif
if set eq 21 then begin dir = '/volume1/scratch/set2/high_T/' & ka_n = '_ka2.24_highT_fe9' & sln = 4 & endif
if set eq 22 then begin dir = '/volume1/scratch/set2/sm/' & ka_n = '_ka2.24_sm_fe9' & sln = 4 & endif
if set eq 23 then begin dir = '/volume1/scratch/set2/193/' & ka_n = '_ka2.24_fe12' & sln = 4 & endif
if set eq 24 then begin dir = '/volume1/scratch/set2/hgres/angs/' & ka_n = '_ka2.24_hgres_fe9' & sln = 5 & endif
;if set eq 24 then begin dir = '/volume1/scratch/set2/hgres/' & ka_n = '_ka2.24_hgres_fe9' & sln = 4 & endif
if set eq 25 then begin dir = '/volume1/scratch/set2/hgres2d/' & ka_n = '_ka2.24_hgres2d.fe9' & sln = 4 & endif
if set eq 26 then begin dir = '/volume1/scratch/set2/hgres2d/' & ka_n = '_ka2.24_hgres2d.fe12' & sln = 4 & endif
if set eq 3 then begin dir='/volume1/scratch/set3/' &  ka_n='_ka1.25_fe9' & sln = 3 & endif
;if set eq 3 then begin dir = '/users/cpa/pantolin/Modeling/cubes/set3/' & ka_n = '_ka1.25_fe9' & sln = 6 & endif

if set ne 25 and set ne 26 then files = file_search(dir+'rslt_slice_*.sav',count=nfiles,/fully_qualify_path) else files = file_search(dir+'rslt_cubes'+ka_n+'*.sav',count=nfiles,/fully_qualify_path)
if set eq 21 or set eq 24 then begin ltn = 24 & lxn = 29 & endif 
if set eq 2 or set eq 23 or set eq 3 then begin ltn = 18 & lxn = 23 & endif
if set eq 22 then begin ltn = 21 & lxn = 26 & endif
if set eq 25 then begin ltn = 30 & lxn = 0 & endif
if set eq 26 then begin ltn = 31 & lxn = 0 & endif

if set eq 25 or set eq 26 then ldimz = 1 else ldimz = 2
if set eq 24 then dimt = 30

restore,files[0]
nwave = n_elements(wave)
;dim00d = (size(image00d_ext))[ldimz]
;if imaging eq 0 then begin
;   imgcube00d = fltarr(dimx,dim00d,nwave,dimt) 
;   endif else begin
;      imgcube00d = fltarr(dimx,dim00d,dimt)
;endelse
if set ne 25 and set ne 26 then dimx = n_elements(gridx) else dimx = 1
if set ne 25 and set ne 26 then begin
;   dim00d = (size(image00d_ext))[ldimz]
   dim30d = (size(image30d_ext))[ldimz]
   dim45d = (size(image45d_ext))[ldimz]
   dim60d = (size(image60d_ext))[ldimz]
   if imaging eq 0 then begin
;      imgcube00d = fltarr(dimx,dim00d,nwave,dimt)
      imgcube30d = fltarr(dimx,dim30d,nwave,dimt)
      imgcube45d = fltarr(dimx,dim45d,nwave,dimt)
      imgcube60d = fltarr(dimx,dim60d,nwave,dimt)
   endif else begin
;      imgcube00d = fltarr(dimx,dim00d,dimt)
      imgcube30d = fltarr(dimx,dim30d,dimt)
      imgcube45d = fltarr(dimx,dim45d,dimt)
      imgcube60d = fltarr(dimx,dim60d,dimt)
   endelse
endif

  for i=0,nfiles-1 do begin
     restore,files[i]
     slice = (strsplit(files[i],'/',/extract))[sln]
     tn = fix(strmid(slice,ltn,3))
     if set ne 25 and set ne 26 then xn = fix(strmid(slice,lxn,3)) 
     if set eq 25 or set eq 26 then xn = 0
     if set eq 24 then xn = fix(strmid(slice,lxn,4))
;     if set ne 24 then imgcube00d[xn,*,*,tn] = image00d_ext else imgcube00d[xn,*,tn] = image00d_ext
     if imaging eq 0 then begin
        if set ne 25 and set ne 26 then begin
;           imgcube00d[xn,*,*,tn] = image00d_ext
           imgcube30d[xn,*,*,tn] = image30d_ext
           imgcube45d[xn,*,*,tn] = image45d_ext
           imgcube60d[xn,*,*,tn] = image60d_ext
           if xn le dimx/2 then begin
              xn2 = dimx-xn-1
;              imgcube00d[xn2,*,*,tn] = imgcube00d[xn,*,*,tn]
              imgcube30d[xn2,*,*,tn] = imgcube30d[xn,*,*,tn]
              imgcube45d[xn2,*,*,tn] = imgcube45d[xn,*,*,tn]
              imgcube60d[xn2,*,*,tn] = imgcube60d[xn,*,*,tn]
           endif
        endif 
     endif else begin
        if set ne 25 and set ne 26 then begin
;           imgcube00d[xn,*,tn] = image00d_ext
           imgcube30d[xn,*,tn] = image30d_ext
           imgcube45d[xn,*,tn] = image45d_ext
           imgcube60d[xn,*,tn] = image60d_ext
           if xn le dimx/2 then begin
              xn2 = dimx-xn-1
;              imgcube00d[xn2,*,tn] = imgcube00d[xn,*,tn]
              imgcube30d[xn2,*,tn] = imgcube30d[xn,*,tn]
              imgcube45d[xn2,*,tn] = imgcube45d[xn,*,tn]
              imgcube60d[xn2,*,tn] = imgcube60d[xn,*,tn]
           endif
        endif
     endelse
;     if set eq 24 then begin
;        if xn le dimx/2 then begin
;           xn2 = dimx-xn-1
;           imgcube00d[xn2,*,tn] = imgcube00d[xn,*,tn]
;        endif
;     endif
     print,string(13b)+' % finished: ',float(i)*100./(nfiles-1),format='(a,f4.0,$)'
  endfor
  if set eq 24 then begin
     for i=0,600 do begin
        for j=0,dimt-1 do begin
           i2 = dimx-i-1
;           imgcube00d[i,*,j] = imgcube00d[601,*,j]
;           imgcube00d[i2,*,j] =imgcube00d[601,*,j]
           imgcube30d[i,*,j] = imgcube30d[601,*,j]
           imgcube30d[i2,*,j] =imgcube30d[601,*,j]
           imgcube45d[i,*,j] = imgcube45d[601,*,j]
           imgcube45d[i2,*,j] =imgcube45d[601,*,j]
           imgcube60d[i,*,j] = imgcube60d[601,*,j]
           imgcube60d[i2,*,j] =imgcube60d[601,*,j]
        endfor
     endfor
  endif
  

  if savcub eq 1 then save,wave,imgcube00d,imgcube30d,imgcube45d,imgcube60d,filename=dir+'imgcubes'+ka_n+'.sav'
end
