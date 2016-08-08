

PRO gridlos, gridx=gridx, gridy=gridy, mua_d=mua_d, velx=velx, vely=vely, dx=dx, dy=dy, n_gridx=n_gridx, n_gridy=n_gridy, ngrid=ngrid, dl=dl, losvel=losvel

if keyword_set(gridx) eq 0 then begin
   print,'gridlos, gridx=gridx, gridy=gridy, mua_d=mua_d, velx=velx, vely=vely, dx=dx, dy=dy, n_gridx=n_gridx, n_gridy=n_gridy, ngrid=ngrid, dl=dl, losvel=losvel'
   return
endif

; Calculates the LOS rays provided a LOS angle, and optionally outputs
; the LOS velocity

; INPUT:
; gridx = (1d float array) grid along x axis
; gridy = (1d float array) grid along y axis
; mua_d = (float) angle between line-of-sight and y-axis
; dx = spatial resolution along gridx-axis
; dy = spatial resolution along gridy-axis

; OUTPUT:
; n_gridx = x-coordinates for points along rays
; n_gridy = y-coordinates for points along rays
; ngrid = provides number and ordering of rays
; dl = distance between rays

; OPTIONAL INPUT:
; velx = (2d float) array of velocity along x
; vely = (2d float) array of velocity along y
; Then the following is ouput: losvel = line-of-sight velocity

; The returned arrays have the values ordered in the following way:
; (n_gridx_<num>[ngrid_<num>[i]:ngrid_<num>[i+1]-1] : x-coordinates for the i-th ray of the <num>-th angle (unscaled)
; (n_gridy_<num>[ngrid_<num>[i]:ngrid_<num>[i+1]-1] : y-coordinates for the i-th ray of the <num>-th angle (unscaled)
; interpol(gridx,indgen(n_elements(gridx)),n_gridx_<num>[ngrid_<num>[i]:ngrid_<num>[i+1]-1]) : x-coordinates on gridx axis for the i-th ray of the <num>-th angle. 
; interpol(gridx,indgen(n_elements(gridx)),n_gridx_<num>[ngrid_<num>[i]:ngrid_<num>[i+1]-1]) : same as above but for the y-coordinates.
; n_elements(ngrid_<num>): number of rays for the <num>-th angle

; mua_d between -90 and 90 deg
if (abs(mua_d) gt 90) then begin
   print,'angles between -90 and 90 degrees'
   return
endif
if (keyword_set(dy) eq 0 or keyword_set(dx) eq 0) then begin dy = 1. & dx = 1. & endif

mua_r=mua_d*!pi/180.
mua_uni = abs(atan(dy/dx*tan(mua_r)))
dimx = float(n_elements(gridx))
dimy = float(n_elements(gridy))
xc = dimx/2.
yc = dimy/2.
if mua_r eq !pi/2 then begin
   yi = yc & ye = yc
   xi = 0 & xe = dimx-1
   x0_p = 0 & x0_m = 0
   xd_p = dimx-1 & xd_m = dimx-1
   i = 1 & ind = 0
   xi_ar = xi
   xe_ar = xe
   yi_ar = yi
   ye_ar = ye
   ls = 1.
   dl = dy
   while(ind eq 0) do begin
      y0_p = yc+float(i)
      yd_p = yc+float(i)
      y0_m = yc-float(i)
      yd_m = yc-float(i)
      if (round(y0_p) gt dimy-1 or round(y0_m) lt 0.) then begin
         ind=1
         xi_ar = [x0_p,xi_ar]
         xe_ar = [xd_p,xe_ar]
         yi_ar = [yi_ar,0.]
         ye_ar = [ye_ar,0.]
      endif else begin
         xi_ar = [x0_p,xi_ar,x0_m]
         xe_ar = [xd_p,xe_ar,xd_m]
         yi_ar = [y0_p,yi_ar,y0_m]
         ye_ar = [yd_p,ye_ar,yd_m]
      endelse
      i = i+1
   endwhile
endif else begin
   dl = sqrt((dx*cos(mua_r))^2+(dy*sin(mua_r))^2)
   ds = sqrt((dx*sin(mua_r))^2+(dy*cos(mua_r))^2)
   if dx lt dy then begin 
;      dl = dx/cos(mua_r)<dy
;      ds = dx/sin(mua_r)<dy
      ls = (ds/dy)<1.
   endif
   if dy lt dx then begin
;      dl = dy/sin(mua_r)<dx
;      ds = dy/cos(mua_r)<dx
      ls = (ds/dx)<1.
   endif
   if dy eq dx then begin
      dl = dx
      ds = dy
      ls = 1.
   endif
   
   xi = (xc-yc*tan(mua_r)*dy/dx)>0.
   xe = (xc+yc*tan(mua_r)*dy/dx)<float(dimx-1)
   if mua_r eq 0. then yi = 0. else yi = (yc-xc*tan(!pi/2.-mua_r)*dx/dy)>0.
   if mua_r eq 0. then ye = float(dimy-1) else ye = (yc+xc*tan(!pi/2.-mua_r)*dx/dy)<float(dimy-1)
   
   xi_ar = xi
   xe_ar = xe
   yi_ar = yi
   ye_ar = ye
   ind = 0
   i = 1
   while(ind eq 0) do begin
      xd_p = (tan(mua_r)*dy/dx*dimy/2.+dimx/2.-float(i)*dl/dx/cos(mua_r))>0.
      xd_m = tan(mua_r)*dy/dx*dimy/2.+dimx/2.+float(i)*dl/dx/cos(mua_r)
      x0_p = -tan(mua_r)*dy/dx*dimy/2.+dimx/2.-float(i)*dl/dx/cos(mua_r)
      x0_m = -tan(mua_r)*dy/dx*dimy/2.+dimx/2.+float(i)*dl/dx/cos(mua_r)
      if mua_r ne 0. then begin
         yd_p = 1./tan(mua_r)*dx/dy*dimx/2.+float(i)*dl/dy/sin(mua_r)+dimy/2.
         yd_m = (1./tan(mua_r)*dx/dy*dimx/2.-float(i)*dl/dy/sin(mua_r)+dimy/2.)>0.
         y0_p = (-1./tan(mua_r)*dx/dy*dimx/2.+float(i)*dl/dy/sin(mua_r)+dimy/2.)>0.
         y0_m = (-1./tan(mua_r)*dx/dy*dimx/2.-float(i)*dl/dy/sin(mua_r)+dimy/2.)>0.
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
;   if i eq 0 then plots,[gridx[xi_ar],gridx[xe_ar]],[gridy[yi_ar],gridy[ye_ar]],color=255
;   plots,[gridx[x0_pp],gridx[xd_pp]],[gridy[y0_pp],gridy[yd_pp]],color=255
;   plots,[gridx[x0_mm],gridx[xd_mm]],[gridy[y0_mm],gridy[yd_mm]],color=255
;   stop
      i = i+1
   endwhile
endelse
;i_sort = sort(xi_ar)
;xi_arr = xi_ar[i_sort]
;xe_arr = xe_ar[i_sort]
;yi_arr = yi_ar[i_sort]
;ye_arr = ye_ar[i_sort]

dimx_n = n_elements(xi_ar)

for i=0., dimx_n-1 do begin
   spline_p, [xi_ar[i],xe_ar[i]], [yi_ar[i],ye_ar[i]], ni_gridx, ni_gridy, interval=ls
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
if n_elements(velx) ne 0 and n_elements(vely) ne 0 then losvel=vely*cos(mua_r)+velx*sin(mua_r)

end
