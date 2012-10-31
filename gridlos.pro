

PRO gridlos, gridx=gridx, gridy=gridy, gridz=gridz, mua_d=mua_d, velx=velx, vely=vely, velz=velz, n_gridx, n_gridy, ngrid, losvel

; INPUT:
; gridx = grid along new x axis (longitudinal direction)
; gridy = grid along new y axis (radial direction)
; gridz = grid along new z axis (radial direction)
; mua_d = angle between line-of-sight and perpendicular to cylinder axis
; velx = (y,z) array of velocity along x
; vely = (y,z) array of velocity along y
; velz = (y,z) array of velocity along z
; OUTPUT:
; n_gridx = grid along new x direction (longitudinal direction)
; n_gridy = grid along new y direction (radial direction)
; ngrid = number of points in depth for each x,y position
; losvel = line-of-sight velocity (in km/s)

; mua_d between -90 and 90 deg
if (abs(mua_d) gt 90) then begin
   print,'angles between -90 and 90 degrees'
   return
endif

mua_r=mua_d*!pi/180.
dimx = n_elements(gridx)
dimy = n_elements(gridy)
dimz = n_elements(gridz)
xc = float(dimx)/2.
yc = float(dimy)/2.

xi = (xc-yc*tan(mua_r))>0.
xe = (xc+yc*tan(mua_r))<float(dimx-1)
if mua_r eq 0. then yi = 0. else yi = (yc-xc*tan(!pi/2.-mua_r))>0.
if mua_r eq 0. then ye = float(dimy-1) else ye = (yc+xc*tan(!pi/2.-mua_r))<float(dimy-1)

xi_ar = xi
xe_ar = xe
yi_ar = yi
ye_ar = ye
ind = 0
i = 1
while(ind eq 0) do begin
   xd_p = (tan(mua_r)*(dimy)/2.+(dimx)/2.-float(i)/cos(mua_r))>0.
   xd_m = tan(mua_r)*(dimy)/2.+(dimx)/2.+float(i)/cos(mua_r)
   x0_p = -tan(mua_r)*(dimy)/2.+(dimx)/2.-float(i)/cos(mua_r)
   x0_m = -tan(mua_r)*(dimy)/2.+(dimx)/2.+float(i)/cos(mua_r)
   if mua_r ne 0. then begin
      yd_p = 1./tan(mua_r)*(dimx)/2.+float(i)/sin(mua_r)+(dimy)/2.
      yd_m = (1./tan(mua_r)*(dimx)/2.-float(i)/sin(mua_r)+(dimy)/2.)>0.
      y0_p = (-1./tan(mua_r)*(dimx)/2.+float(i)/sin(mua_r)+(dimy)/2.)>0.
      y0_m = (-1./tan(mua_r)*(dimx)/2.-float(i)/sin(mua_r)+(dimy)/2.)>0.
   endif else begin
      yd_p = float(dimy-1)
      yd_m = float(dimy-1)
      y0_p = 0.
      y0_m = 0.
   endelse
   if x0_p lt 0. then begin
      x0_pp = 0. & y0_pp = y0_p
   endif else begin
      x0_pp = x0_p & y0_pp = 0.
   endelse
   if x0_m lt 0. then begin      
      x0_mm = 0. & y0_mm = y0_m
   endif else begin
      x0_mm = x0_m & y0_mm = 0.
   endelse
   if xd_p gt float(dimx-1) then begin
      xd_pp = float(dimx-1) & yd_pp = yd_p
   endif else begin
      xd_pp = xd_p & yd_pp = float(dimy-1)
   endelse
   if xd_m gt float(dimx-1) then begin
      xd_mm = float(dimx-1) & yd_mm = yd_m
   endif else begin
      xd_mm = xd_m & yd_mm = float(dimy-1)
   endelse      
   if (round(x0_pp) le 0. and round(xd_pp) le 0.) then begin
      ind=1
      if round(x0_pp) eq 0. and round(x0_mm) eq dimx then begin
         xi_ar = [x0_pp,xi_ar]
         xe_ar = [xd_pp,xe_ar]
         yi_ar = [y0_pp,yi_ar]
         ye_ar = [yd_pp,ye_ar]
      endif      
   endif else begin
      xi_ar = [x0_pp,xi_ar,x0_mm]
      xe_ar = [xd_pp,xe_ar,xd_mm]
      yi_ar = [y0_pp,yi_ar,y0_mm]
      ye_ar = [yd_pp,ye_ar,yd_mm]
   endelse
   i = i+1
endwhile

;i_sort = sort(xi_ar)
;xi_arr = xi_ar[i_sort]
;xe_arr = xe_ar[i_sort]
;yi_arr = yi_ar[i_sort]
;ye_arr = ye_ar[i_sort]

dimx_n = n_elements(xi_ar)
;stop
for i=0, dimx_n-1 do begin
   spline_p, [xi_ar[i],xe_ar[i]], [yi_ar[i],ye_ar[i]], ni_gridx, ni_gridy, interval=1
   nigrid = n_elements(ni_gridx)
   if i eq 0 then begin
      n_gridx = ni_gridx
      n_gridy = ni_gridy
      ngrid = [0,long(nigrid)]
   endif else begin
      n_gridx = [n_gridx,ni_gridx]
      n_gridy = [n_gridy,ni_gridy]
      ngrid = [ngrid, long(ngrid[i]+nigrid)]
   endelse
endfor

; calculate line-of-sight velocity array:
losvel=vely*cos(mua_r)+velx*sin(mua_r)


; construct triangulation of the given points in a plane:
;triangulate, gridx, gridy, triangles, boundary

end
