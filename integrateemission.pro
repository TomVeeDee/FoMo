pro integrateemission,emis=emis,n_gridx=n_gridx,n_gridy=n_gridy,ngrid=ngrid,wave=wave,w0=w0,direction=direction,losvel=losvel,image

; INPUT:	
; emis = output of lineongrid.pro , array of emission (x,y,(z),lambda)
; wave = wavelength array from lineongrid.pro
; n_gridx = output of gridlos.pro, grid along new x direction
; (longitudinal direction: old z axis)
; n_gridy = output of gridlos.pro, grid along new y direction 
; (radial direction)
; ngrid = output of gridlos.pro, number of points in depth for each
; x,y position
; direction = direction of the integration (1 = x, 2 = y, 3 = z, 4 = mua_d angle)
; losvel = output of gridlos.pro, line-of-sight velocity (in km/s)
; OUTPUT:
; image	= intensity along line-of-sight (new x, (new y),lambda) = (nz,ny,nwave)

sizes=size(emis)
dims=sizes[0]-1
nx=sizes[1]
ny=sizes[2]
if dims eq 2 then nz=1 else nz = sizes[3]
if dims eq 2 then nwave = sizes[3] else nwave=sizes[4]

c=299792000./1.e5
emission = emis
doppleremission=emission

; calculating doppler shifts point by point:
;for i=0,nx-1 do begin
;	print,'doing i=',i
;	for j=0,ny-1 do begin
;		if (dims eq 2) then begin
;			line=reform(doppleremission[i,j,*])
;			newwave=wave+losvel[i,j]/c*mean(wave)
;			newline=interpol(line,wave,newwave,/spline)
;			doppleremission[i,j,*]=newline
;		endif
;		if (dims eq 3) then begin
;			for k=0,nz-1 do begin
;				line=reform(doppleremission[i,j,k,*])
;				newwave=wave+losvel[i,j,k]/c*mean(wave)
;				newline=interpol(line,wave,newwave,/spline)
;				doppleremission[i,j,k,*]=newline
;			endfor
;		endif
;	endfor
;endfor

; calculate doppler shifts through binning of velocity matrix
bsize = 0.005
nhlosvel = histogram(losvel,binsize=bsize,locations=histvel,reverse_indices=R)
nhist = n_elements(nhlosvel)
nlosvel = fltarr(nx,ny,nz,nwave)
for i=0,nwave-1 do nlosvel[*,*,*,i]=losvel
nnhlosvel = histogram(nlosvel,binsize=bsize,reverse_indices=Rn)
if dims eq 2 then lemmx = ([max(emission[nx/2,ny/2,*]),!c])[1]
if dims eq 3 then lemmx = ([max(emission[nx/2,ny/2,nz/2,*]),!c])[1]
emipk = reform(emission[0,*,*,lemmx])

for i=0,nhist-1 do begin
   losv = histvel[i]
   newwave = wave+losv/c*w0;mean(wave)
   if R[i] ne R[i+1] then begin
      nR = n_elements(R[R[i]:R[i+1]-1])
      emi = emipk[R[R[i]:R[i+1]-1]]
      indx = array_indices(emipk,R[R[i]:R[i+1]-1])
      numemi = min([n_elements(uniq(emi)),10000])
      hemi = histogram(emi,nbins=numemi,locations=lhemi,reverse_indices=Re)
      nhemi = n_elements(hemi)
      for j=0,nhemi-1 do begin
         if Re[j] ne Re[j+1] then begin
            dopemi=interpol(emission[0,indx[0,Re[Re[j]]],indx[1,Re[Re[j]]],*],wave,newwave,/spline)
            for k=0,hemi[j]-1 do doppleremission[0,indx[0,Re[Re[j]+k]],indx[1,Re[Re[j]+k]],*]=dopemi
         endif
      endfor
   endif
endfor


; then sum up over one direction.

if (direction le dims) then begin
	image=total(doppleremission,direction)
endif else begin
   if direction eq 4 then begin
      if dims eq 2 then image = dblarr(n_elements(ngrid)-1,nwave) else image = dblarr(nx,n_elements(ngrid)-1,nwave)
      for i=0, n_elements(ngrid)-2 do begin
         if (dims eq 2) then begin
            for j=0,nwave-1 do image[i,j] = total(interpolate(reform(doppleremission[*,*,j]),n_gridx[ngrid[i]:ngrid[i+1]-1],n_gridy[ngrid[i]:ngrid[i+1]-1]))
         endif
         if (dims eq 3) then begin
            for k=0,nx-1 do begin
               for j=0,nwave-1 do image[k,i,j] = total(interpolate(reform(doppleremission[k,*,*,j]),n_gridy[ngrid[i]:ngrid[i+1]-1],n_gridx[ngrid[i]:ngrid[i+1]-1]),/double)
            endfor
         endif
         print,string(13b)+' % finished: ',float(i)*100./(n_elements(ngrid)-2),format='(a,f4.0,$)'
      endfor
	print,' '
   endif else begin
      print,'direction for LOS integration higher than dimension'
      image=doppleremission
   endelse

endelse

end
