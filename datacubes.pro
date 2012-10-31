
pro datacubes,gridx,gridy,dimz,dimw,dimt,vr_t,vz_t,rtot_t,te_t,vr_cube,vz_cube,te_cube,rh_cube

if n_params(0) lt 1 then begin
   print,'datacubes,gridx,gridy,dimz,dimw,dimt,vr_t,vz_t,rtot_t,te_t,vr_cube,vz_cube,te_cube,rh_cube'
   return
endif

sav = 0
r0 = gridx[-1]/2.
dimx = n_elements(gridx)
dimy = n_elements(gridy)
dimt0 = 1 & jj = 0
dimw0 = 1 & ii = 14
vr_cube=fltarr(dimw0,dimx,dimy,dimz,dimt0)
vz_cube=fltarr(dimw0,dimx,dimy,dimz,dimt0)
te_cube=fltarr(dimw0,dimx,dimy,dimz,dimt0)
rh_cube=fltarr(dimw0,dimx,dimy,dimz,dimt0)
dist_cube_b=make_array(dimw0,dimx,dimy,dimz,dimt0)
dist_cube=fltarr(dimx,dimy,dimz)

; (frequency, x, y, z, time)

for i=0,dimx-1 do for j=0,dimy-1 do dist_cube[i,j,*] = sqrt((r0-gridx[i])^2+(r0-gridy[j])^2)

for i=0,dimw0-1 do for j=0,dimt0-1 do dist_cube_b[i,*,*,*,j] = dist_cube
mxdist = sqrt((gridx[-1]-r0)^2+(gridy[-1]-r0)^2)

for j=0,dimt-1 do begin
   for l=0,round(mxdist) do begin
      cyl = where(round(dist_cube_b) eq l)
      dis = where(round(reform(dist_cube_b[0,*,*,0,0])) eq l,ndist)
      col3dvr = fltarr(ndist*dimz)
      col3dvz = fltarr(ndist*dimz)
      col3dte = fltarr(ndist*dimz)
      col3drh = fltarr(ndist*dimz)
      for i=0,dimw0-1 do begin
         if l lt max(abs(dimx-r0)) then colvr = reform(vr_t[ii,r0-l,*,j]) else colvr = reform(vr_t[ii,0,*,j])
         for k=0,dimz-1 do col3dvr[k*ndist:(k+1)*ndist-1]=replicate(colvr[k],ndist)
         if l lt max(abs(dimx-r0)) then colvz = reform(vz_t[ii,r0-l,*,j]) else colvz = reform(vz_t[ii,0,*,j])
         for k=0,dimz-1 do col3dvz[k*ndist:(k+1)*ndist-1]=replicate(colvz[k],ndist)
         if l lt max(abs(dimx-r0)) then colte = reform(te_t[ii,r0-l,*,j]) else colte = reform(te_t[ii,0,*,j])
         for k=0,dimz-1 do col3dte[k*ndist:(k+1)*ndist-1]=replicate(colte[k],ndist)
         if l lt max(abs(dimx-r0)) then colrh = reform(rtot_t[ii,r0-l,*,j]) else colrh = reform(rtot_t[ii,0,*,j])
         for k=0,dimz-1 do col3drh[k*ndist:(k+1)*ndist-1]=replicate(colrh[k],ndist)

         ; if dimt0=1 then set jj=1 else if dimt0 = dimt then jj=j

         vr_cube[cyl[i+jj*n_elements(cyl)/dimt0:n_elements(cyl)/dimt0+jj*n_elements(cyl)/dimt0-1:dimw0]] = col3dvr
         vz_cube[cyl[i+jj*n_elements(cyl)/dimt0:n_elements(cyl)/dimt0+jj*n_elements(cyl)/dimt0-1:dimw0]] = col3dvz
         te_cube[cyl[i+jj*n_elements(cyl)/dimt0:n_elements(cyl)/dimt0+jj*n_elements(cyl)/dimt0-1:dimw0]] = col3dte
         rh_cube[cyl[i+jj*n_elements(cyl)/dimt0:n_elements(cyl)/dimt0+jj*n_elements(cyl)/dimt0-1:dimw0]] = col3drh
      endfor
   endfor
stop
   if sav eq 1 then begin
      vr_cube = reform(vr_cube)
      vz_cube = reform(vz_cube)
      te_cube = reform(te_cube)
      rh_cube = reform(rh_cube)
      save,vr_cube,vz_cube,te_cube,rh_cube,filename='cubes_'+string(j,format="(i3.3)")+'.sav'
   endif
endfor

;      for k=0,dimx-1 do begin
;         for l=0,dimy-1 do begin
;            if (gridx[k]+gridy[l]) eq 0 then vr_cube[i,k,l,*,j]=vr_t[i,k,*,j] else vr_cube[i,k,l,*,j]=(abs(gridx[k]-r;0)*vr_t[i,k,*,j]+abs(gridy[l]-r0)*vr_t[i,l,*,j])/(gridx[k]+gridy[l])
;         endfor
;      endfor
;   endfor
;endfor

end


