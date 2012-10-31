
pro mkimgcubes, set=set, gridx=gridx, dimt=dimt, dim00d=dim00d, dim30d=dim30d, dim45d=dim45d, dim60d=dim60d, gridy, gridz_ext, wave, imgcube00d,imgcube30d,imgcube45d,imgcube60d

; Produces cubes from sav files produced by prl_slices.pro
; INPUT:
; set = defines model:
;     = 2: 'base model - 171' corresponds to ka = 2.24, line Fe IX 171
;     = 21: 'high T model' base model with high external temperature
;     = 22: 'smooth model' base model with smooth density profile across cylinder
;     = 23: 'base model - 193' base model with line Fe XII 193
;     = 3 : 'long lambda model' corresponds to ka = 1.25
;     = 24: 'high_res' corresponds to base model with high resolution
; gridx = grid along x axis (radial direction)
; dimt = number of points in time dimension
; dim00d = number of points in longitudinal direction in new grid with angle=00
; dim30d = number of points in longitudinal direction in new grid with angle=30
; dim45d = number of points in longitudinal direction in new grid with angle=45
; dim60d = number of points in longitudinal direction in new grid with angle=60

  savcub = 1
  dimx = n_elements(gridx) & nwave = 100
;  if set eq 2 then begin
;     dim00d = 342 & dim30d = 398 & dim45d = 386 & dim60d = 348
;  endif
;  if set eq 3 or set eq 31 then begin
;     dim00d = 404 & dim30d = 452 & dim45d = 430 & dim60d = 379
;  endif
  imgcube00d = fltarr(dimx,dim00d,nwave,dimt)
  imgcube30d = fltarr(dimx,dim30d,nwave,dimt)
  imgcube45d = fltarr(dimx,dim45d,nwave,dimt)
  imgcube60d = fltarr(dimx,dim60d,nwave,dimt)

if set eq 2 then begin dir = '/volume1/scratch/set2/' & ka_n = '_ka2.24_fe9' & sln = 3 & endif
;if set eq 2 then begin dir = '/users/cpa/pantolin/Modeling/cubes/set2/' & ka_n = '_ka2.24_fe9' & sln = 6 & endif
;if set eq 21 then begin dir = '/users/cpa/pantolin/Modeling/cubes/set2/highT/' & ka_n = '_ka2.24_highT_fe9' & sln = 7 & endif
if set eq 21 then begin dir = '/volume1/scratch/set2/high_T/' & ka_n = '_ka2.24_highT_fe9' & sln = 4 & endif
if set eq 22 then begin dir = '/volume1/scratch/set2/sm/' & ka_n = '_ka2.24_sm_fe9' & sln = 4 & endif
if set eq 23 then begin dir = '/volume1/scratch/set2/193/' & ka_n = '_ka2.24_fe12' & sln = 4 & endif
if set eq 24 then begin dir = '/volume1/scratch/set2/hgres/' & ka_n = '_ka2.24_hgres_fe9' & sln = 4 & endif
  if set eq 3 then begin dir='/volume1/scratch/set3/' &  ka_n='_ka1.25_fe9' & sln = 3 & endif
;if set eq 3 then begin dir = '/users/cpa/pantolin/Modeling/cubes/set3/' & ka_n = '_ka1.25_fe9' & sln = 6 & endif

   files = file_search(dir+'rslt_slice_*.sav',count=nfiles,/fully_qualify_path)
if set eq 21 or set eq 24 then begin ltn = 24 & lxn = 29 & endif 
if set eq 2 or set eq 23 or set eq 3 then begin ltn = 18 & lxn = 23 & endif
if set eq 22 then begin ltn = 21 & lxn = 26 & endif

  for i=0,nfiles-1 do begin
     restore,files[i]
     slice = (strsplit(files[i],'/',/extract))[sln]
     tn = fix(strmid(slice,ltn,3))
     xn = fix(strmid(slice,lxn,3))
     imgcube00d[xn,*,*,tn] = image00d_ext
     imgcube30d[xn,*,*,tn] = image30d_ext
     imgcube45d[xn,*,*,tn] = image45d_ext
     imgcube60d[xn,*,*,tn] = image60d_ext
     if xn le dimx/2 then begin
        xn2 = dimx-xn-1
        imgcube00d[xn2,*,*,tn] = imgcube00d[xn,*,*,tn]
        imgcube30d[xn2,*,*,tn] = imgcube30d[xn,*,*,tn]
        imgcube45d[xn2,*,*,tn] = imgcube45d[xn,*,*,tn]
        imgcube60d[xn2,*,*,tn] = imgcube60d[xn,*,*,tn]
     endif
     print,string(13b)+' % finished: ',float(i)*100./(nfiles-1),format='(a,f4.0,$)'
  endfor
  print, ' '
  if savcub eq 1 then save,wave,imgcube00d,imgcube30d,imgcube45d,imgcube60d,filename=dir+'imgcubes'+ka_n+'.sav'
end
