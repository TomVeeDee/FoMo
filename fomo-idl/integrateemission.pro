pro integrateemission,emission=emission,logt=logt,n_gridx=n_gridx,n_gridy=n_gridy,ngrid=ngrid,wave=wave,w0=w0,direction=direction,losvel=losvel,imaging=imaging,watom=watom,image=image,dl=dl,ds=ds,d_perp=d_perp,wayemi=wayemi,channel=channel,silent=silent

; Calculates intensity by integrating emissivity along a given
; line-of-sight

; INPUT:	
; emission = output of lineongrid_goft_tab.pro: emission_goft, 2d array of emissivity
; wave = wavelength array from lineongrid_goft_tab.pro
; n_gridx, n_gridy, ngrid = output of gridlos.pro
; direction = (int) direction of the integration (1 = x, 2 = y, 3 = z, 4 = mua_d angle)
; losvel = (2d array) output of gridlos.pro, line-of-sight velocity
; (in km/s unit). Only necessary for spectral calculations (no imaging)
; dl = distance between rays
; ds = distance along rays
; d_perp = distance between fomo slices (resolution in 3rd axis: if
; fomo slice is in (x,y) plane then d_perp = dz, where dz is
; resolution along z axis)
  
; OPTIONAL:
; set keyword imaging or channel (-> wayemi = 5) for no doppler shift
; calculation (imaging case) (default = 1).  

; OUTPUT:
; image	= (2d float array) intensity along line-of-sight (n_elements(ngrid),nwave)

; solid angle of 1 simulation pixel in arcsec

if keyword_set(channel) then apix = (dl/1.e5/715.)*(d_perp/1.e5/715.)*(!pi/180./3600.)^2 else apix = 1

sizes=size(emission)
dims = sizes[0]
nx = sizes[1]
ny = sizes[2]
if dims eq 3 then nz = sizes[3] else nz = 1
if ~keyword_set(imaging) and ~keyword_set(channel) then begin
   nwave = n_elements(wave)
endif else begin
   nwave = 1
endelse
proton=1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
c=299792000.d

if dims eq 3 then doppleremission = fltarr(nx,ny,nz,nwave) else doppleremission = fltarr(nx,ny,nwave)

if wayemi ne 5 and ~keyword_set(channel) and ~keyword_set(imaging) then begin

; calculate doppler shifts through binning of velocity matrix
   bsize = 0.5 ; bins of 0.5 km/s
   nhlosvel = histogram(losvel,binsize=bsize,locations=histvel,reverse_indices=R)
   nhist = n_elements(nhlosvel)

   if wayemi ne 4 then begin
      if dims eq 3 then begin
         lemmx = ([max(emission[nx/2,ny/2,nz/2,*]),!c])[1]
         if lemmx eq 0. then lemmx=nwave/2.
         emipk = reform(emission[0,*,*,lemmx])
      endif else begin
         lemmx = ([max(emission[nx/2,ny/2,*]),!c])[1]
         if lemmx eq 0. then lemmx=nwave/2.
         emipk = reform(emission[*,*,lemmx])
      endelse
   endif else begin
      emipk = emission
   endelse

   for i=0.,nhist-1 do begin
      losv = histvel[i]
      newwave = wave+losv*1.e3/c*w0  ;mean(wave)
      wdop = w0+double(losv*1.e3/c*w0)
      if R[i] ne R[i+1] then begin
         nR = n_elements(R[R[i]:R[i+1]-1])
         emi = emipk[R[R[i]:R[i+1]-1]]
         indx = array_indices(emipk,R[R[i]:R[i+1]-1])
         numemi = min([n_elements(uniq(emi)),10000])
         hemi = histogram(emi,nbins=numemi,locations=lhemi,reverse_indices=Re)
         nhemi = n_elements(hemi)
         for j=0.,nhemi-1 do begin
            if (Re[j] ne Re[j+1] and lhemi[j] ne 0.) then begin
               mnlgt = mean(logt[indx[0,Re[Re[j]:Re[j+1]-1]],indx[1,Re[Re[j]:Re[j+1]-1]]])
               sigma = sqrt(kboltz/proton/watom*w0^2/c^2*10^(mnlgt))
               prms =[lhemi[j],wdop,sigma]
               dopemi = gaussian(wave,prms)/sigma/sqrt(2*!Pi)
               if dims eq 3 then begin
;                  dopemi=interpol(emission[0,indx[0,Re[Re[j]]],indx[1,Re[Re[j]]],*],wave,newwave,/spline) 
                  for k=0.,hemi[j]-1 do doppleremission[0,indx[0,Re[Re[j]+k]],indx[1,Re[Re[j]+k]],*]=dopemi
               endif else begin 
                  
;                  dopemi=interpol(emission[indx[0,Re[Re[j]]],indx[1,Re[Re[j]]],*],wave,newwave,/spline)
                  for k=0.,hemi[j]-1 do doppleremission[indx[0,Re[Re[j]+k]],indx[1,Re[Re[j]+k]],*]=dopemi
               endelse
            endif
         endfor
      endif
   endfor
endif else begin
   doppleremission = emission
endelse

; then sum up over rays:

if (direction le dims) then begin
   image=total(doppleremission,direction)*ds*apix
   if direction eq 1 then image = reverse(image,1)
endif else begin
   if direction eq 4 then begin
      if dims eq 2 then image = dblarr(n_elements(ngrid)-1,nwave) else image = dblarr(nx,n_elements(ngrid)-1,nwave)
      for i=0., n_elements(ngrid)-2 do begin
         if (dims eq 2) then begin
            for j=0,nwave-1 do image[i,j] = total(interpolate(reform(doppleremission[*,*,j]),n_gridx[ngrid[i]:ngrid[i+1]-1],n_gridy[ngrid[i]:ngrid[i+1]-1]))*ds*apix
         endif
         if (dims eq 3) then begin
            for k=0,nx-1 do begin
               for j=0,nwave-1 do image[k,i,j] = total(interpolate(reform(doppleremission[k,*,*,j]),n_gridx[ngrid[i]:ngrid[i+1]-1],n_gridy[ngrid[i]:ngrid[i+1]-1]),/double)*ds*apix
            endfor
         endif
         if ~keyword_set(silent) then print,string(13b)+' % finished: ',float(i)*100./(n_elements(ngrid)-2),format='(a,f4.0,$)'
      endfor
      if ~keyword_set(silent) then print,' '
   endif else begin
      print,'direction for LOS integration higher than dimension'
   endelse
endelse

end
