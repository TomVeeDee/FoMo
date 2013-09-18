
pro rays_pick, x0=x0, y0=y0, n_gridx=n_gridx, n_gridy=n_gridy, ngrid=ngrid, mua=mua, nray

if keyword_set(x0) eq 0 then begin
   print,'x0=x0, y0=y0, n_gridx=n_gridx, n_gridy=n_gridy, ngrid=ngrid, mua=mua, nray'
   return
endif

; Returns the ray number crossing at coordinates (x0,y0) 

; INPUT:
; x0: (float) x-coordinate crossing where ray is desired
; y0: (float) y-coordinate crossing where ray is desired
; n_gridx, n_gridy, ngrid: defined in gridlos.pro
; mua: line-of-sight angle

; OUTPUT:
; nray: ray number crossing at coordinates (x0,y0)

; Eaxmple of use once nray is obtained:
; (n_gridx[ngrid[nray]:ngrid[nray+1]-1] : x-coordinates for nray 
; (n_gridy[ngrid[nray]:ngrid[nray+1]-1] : y-coordinates for nray 

  ang = string(mua,format="(i2.2)")
  locxy = fltarr(n_elements(ngrid)-1,2)
  for i=0.,n_elements(ngrid)-2 do begin
     locxy[i,*]=[min(sqrt((n_gridy[ngrid[i]:ngrid[i+1]-1]-y0)^2+(n_gridx[ngrid[i]:ngrid[i+1]-1]-x0)^2)),!c]
  endfor
  nray = ([min(locxy[*,0]),!c])[1]
  print,'LOS '+ang+' ray passing through (x0,y0)',nray
end
