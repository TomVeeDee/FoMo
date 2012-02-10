pro integrateemission,emission,rho,grid,wave,direction,losvel,image

; input:	emission	output of lineongrid
;		rho		density
;		wave		wavelengths from lineongrid
;		grid		the grid of the emission data
;		direction	direction of the integration (1 = x, 2 = y, 3 = z)
;		losvel		line-of-sight velocity (in m/s), has the same dimensions as emission[*,*,*,0]
; output:	image		3D cube containing a spectroscopic image, dimensions e.g. ny x nz x nwave

; first do all the doppler shifts
; I'll use the naive method: shift wavelength, interpolate the original data to the new wavelength, and load the data into the emission matrix
; this has to use a long for-loop though. :-(

sizes=size(emission)
dims=sizes[0]-1
nx=sizes[1]
ny=sizes[2]
nz=sizes[3]
if (dims eq 2) then nz=1
nwave=sizes[4]

c=299792000.

doppleremission=emission

for i=0,nx-1 do begin
	print,'doing i=',i
	for j=0,ny-1 do begin
		if (dims eq 2) then begin
			line=reform(doppleremission[i,j,*])
			newwave=wave+losvel[i,j]/c*mean(wave)
			newline=interpol(line,wave,newwave,/spline)
			doppleremission[i,j,*]=newline
		endif
		if (dims eq 3) then begin
			for k=0,nz-1 do begin
				line=reform(doppleremission[i,j,k,*])
				newwave=wave+losvel[i,j,k]/c*mean(wave)
				newline=interpol(line,wave,newwave,/spline)
				doppleremission[i,j,k,*]=newline
			endfor
		endif
	endfor
endfor

for i=0,nwave-1 do begin
	if (dims eq 2) then doppleremission[*,*,i]=doppleremission[*,*,i]*rho^2
	if (dims eq 3) then doppleremission[*,*,*,i]=doppleremission[*,*,*,i]*rho^2
endfor

; then sum up over one direction.

if (direction le dims) then begin
	image=total(doppleremission,direction)
endif else begin
	print,'direction for LOS integration higher than dimension'
	image=doppleremission
endelse

end
