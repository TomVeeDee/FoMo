
pro degrad_res,oslice=oslice,wave=wave,w0=w0,dlx=dlx,dly=dly,gridx=gridx,gridy=gridy,dt=dt,time=time,xres=xres,fwhm_xres=fwhm_xres,yres=yres,fwhm_yres=fwhm_yres,vres=vres,cad=cad,noise=noise,effarea=effarea,dslice_cgn=dslice_cgn,grx_cg=grx_cg,gry_cg=gry_cg,wav_cg=wav_cg,imaging=imaging,spdata=spdata,wfluc=wfluc

; INPUT:

; oslice : data in format: (x,wavelength,time) or
;                          (x,y,wavelength,time) or
;                          (x,y,time) or
;                          (x,time)
;         in units of erg cm^-2 s^-1 sr^-1
  
; x(y)res, x(y)res_fwhm : in units of arcsec
; dlx : x-grid size in arcsec
; (dly : y-grid size in arcsec)
; cad : cadence of target instrument
; dt = time unit of data in sec
; time = time array in sec
; effarea : effective area of target instrument in specific emission
;           line (in cm^2)
; For imaging output only set keyword 'imaging' and
;      and provide centre wavelength of emission w0 (Angstrom)

; If provided data has spectral data then set keyword 'spdata'
; For spectral data (no /spdata):
;     w0 : centre wavelength of emission line (in Angstrom)
;     if no /imaging then vres : in units of km/s
;     wave : wavelength array of emission line
  
; OPTIONAL:
; fwhm_xres : target x-resolution of the PSF (taken as the FWHM)
;               default: fwhm_xres = x_res*2
; fwhm_yres : target y-resolution of the PSF (taken as the FWHM)
;               default: fwhm_yres = y_res*2
; gridx : x-grid
; gridy : y-grid
; wfluc : amplitude of random fluctuations for noise at each
;        wavelength position. Default is 0.1 (=10% of intensity)

; OUTPUT:
;
; dslice_cgn : degraded data, with dimensions
;     oslice                  ->      dslice_cgn
;    if (dimx,dimy,dimw,dimt) -> (nx_cg,ny_cg,nw_cg,num_t))
; if gridx is present then grx_cg is the new degraded x-grid
; if gridy is present then gry_cg is the new degraded y-grid
  
  h_cgs = 6.6261e-27
  c_cgs = 2.9979e10

  dims = size(oslice)
  if dims[0] eq 4 and ~keyword_set(imaging) then begin
; old_array(dimx,dimy,dimw,dimt) -> new_array(new_dimx,new_dimy,new_dimw,new_dimt)
     dimx = dims[1]
     dimy = dims[2]
     dimw = dims[3]
     dimt = dims[4]
     locw = 3
     slice = oslice
  endif
  if dims[0] eq 4 and keyword_set(imaging) then begin
; old_array(dimx,dimy,dimw,dimt) -> new_array(new_dimx,new_dimy,new_dimt)
     dimx = dims[1]
     dimy = dims[2]
     dimw = 0
     dimt = dims[4]
     locw = 3
     slice = total(oslice,locw)
  endif
  if dims[0] eq 3 and ~keyword_set(imaging) and keyword_set(spdata) then begin
; old_array(dimx,dimw,dimt) -> new_array(new_dimx,new_dimw,new_dimt)
     dimx = dims[1]
     dimy = 0
     dimw = dims[2]
     dimt = dims[3]
     locw = 2
     slice = oslice
  endif
  if dims[0] eq 3 and ~keyword_set(spdata) then begin
; old_array(dimx,dimy,dimt) -> new_array(new_dimx,new_dimy,new_dimt)
     dimx = dims[1]
     dimy = dims[2]
     dimw = 0
     dimt = dims[3]
     slice = oslice
  endif
  if dims[0] eq 3 and keyword_set(imaging) and keyword_set(spdata) then begin
; old_array(dimx,dimw,dimt) -> new_array(new_dimx,new_dimt)
     dimx = dims[1]
     dimy = 0
     dimw = 0
     dimt = dims[3]
     locw = 2
     slice = total(oslice,locw)
  endif
  if dims[0] eq 2 then begin
; old_array(dimx,dimt) -> new_array(new_dimx,new_dimt)
     dimx = dims[1]
     dimy = 0
     dimw = 0
     dimt = dims[2]
     slice = oslice
  endif
  
  c = 299792.458
  if dimw ne 0 then begin
     cf = h_cgs*c_cgs/(wave*1.e-8)
     dw = mean(wave[1:*]-wave[0:*])
     nwave = n_elements(wave)
     wres = vres*w0/c
     numpix_w = wres/dw
     nw_cg = round(nwave/numpix_w)
     wav_cg = congrid(wave,nw_cg,/center,/interp)
     dw_cg = mean(wav_cg[1:*]-wav_cg[0:*])
  endif else begin
     if n_elements(w0) ne 0 then cf = h_cgs*c_cgs/(w0*1.e-8)
     if dims[0] eq 3 then dw = mean(wave[1:*]-wave[0:*]) else dw = 1
  endelse

  if ~keyword_set(wfluc) and dimw ne 0 then wfluc = 0.1

  apix = xres^2*2.35d-11

  if n_elements(vres) ne 0 then begin
  endif

  if ~keyword_set(fwhm_xres) then fwhm_xres = xres*2
  if ~keyword_set(fwhm_yres) and n_elements(yres) ne 0. then fwhm_yres = yres*2
  numpix_x = xres/dlx
  numpix_x_fwhm = fwhm_xres/dlx
  nx_cg = round(dimx/numpix_x)
  if n_elements(gridx) ne 0 then grx_cg = congrid(gridx,nx_cg,/center,/interp)

  if dimy ne 0 then begin
     numpix_y = yres/dly
     numpix_y_fwhm = fwhm_yres/dly
     ny_cg = round(dimy/numpix_y)     
     if n_elements(gridy) ne 0 then gry_cg = congrid(gridy,ny_cg,/center,/interp)
  endif

  if ~keyword_set(cad) then begin
     num_t = dimt
     cad_eff = dt
     time_cad = time
  endif else begin
     num_t = round(dimt/(cad/dt))
     time_cad = congrid(time,num_t,/interp)
     cad_eff = time_cad[1]
  endelse

  if dims[0] eq 4 and ~keyword_set(imaging) then begin
     dslice_cg = fltarr(nx_cg,ny_cg,nw_cg,num_t)
  endif
  if dims[0] eq 4 and keyword_set(imaging) then begin
     dslice_cg = fltarr(nx_cg,ny_cg,num_t)
  endif
  if dims[0] eq 3 and ~keyword_set(imaging) then begin
     dslice_cg = fltarr(nx_cg,nw_cg,num_t)
  endif
  if dims[0] eq 3 and keyword_set(imaging) and keyword_set(spdata) then begin
     dslice_cg = fltarr(nx_cg,num_t)
  endif
  if dims[0] eq 3 and ~keyword_set(spdata) then begin
     dslice_cg = fltarr(nx_cg,ny_cg,num_t)
  endif
  if dims[0] eq 3 and ~keyword_set(imaging) and keyword_set(spdata) then begin
     dslice_cg = fltarr(nx_cg,nw_cg,num_t)
  endif
  if dims[0] eq 2 then begin
     dslice_cg = fltarr(nx_cg,num_t)
  endif
     
  if n_elements(cad) ne 0 then begin
     if dimw ne 0 then begin
        if dimy eq 0 then begin
           slice_t = fltarr(dimx,dimw,num_t)
           for i=0,dimx-1 do slice_t[i,*,*] = frebin(reform(slice[i,*,*]),nwave,num_t,/total)/cad_eff
        endif else begin
           slice_t = fltarr(dimx,dimy,dimw,num_t)
           for j=0,dimy-1 do begin
              for i=0,dimx-1 do begin
                 slice_t[i,j,*,*] = frebin(reform(slice[i,j,*,*]),nwave,num_t,/total)/cad_eff
              endfor
              print,string(13b)+' % finished: ',float(j)*100./((dimy-1)>1),format='(a,f4.0,$)'
           endfor
        endelse
     endif else begin
        if dimy eq 0 then begin
           slice_t = frebin(slice,dimx,num_t,/total)/cad_eff
        endif else begin
           slice_t = fltarr(dimx,dimy,num_t)
           for i=0,dimx-1 do slice_t[i,*,*] = frebin(reform(slice[i,*,*]),dimy,num_t,/total)/cad_eff
        endelse
     endelse        
  endif else begin
     slice_t = slice/cad_eff
  endelse

  if dimw ne 0 then begin
     if dimy ne 0 then begin
        slice_s = slice_t*0.
        slice_w = slice_t*0.
        dslice_cgw = fltarr(nx_cg,ny_cg,dimw,num_t)
        for j=0,dimw-1 do begin
           for i=0,num_t-1 do begin
              if max(slice_t[*,*,j,i]) gt 0 then begin
                 slice_s[*,*,j,i] = filter_image(reform(slice_t[*,*,j,i]),FWHM=[numpix_x_fwhm,numpix_y_fwhm],/all)
              endif
           endfor
           print,string(13b)+' % finished: ',float(j)*100./((dimw-1)>1),format='(a,f4.0,$)'
        endfor
        for i=0,dimx-1 do begin
           for j=0,dimy-1 do begin
              if max(slice_s[i,j,*,*]) gt 0 then begin
                 slice_w[i,j,*,*] = filter_image(reform(slice_s[i,j,*,*]),FWHM=[numpix_w,1],/all)
              endif
           endfor
           print,string(13b)+' % finished: ',float(i)*100./((dimx-1)>1),format='(a,f4.0,$)'
        endfor
        for i=0,num_t-1 do begin
           for j=0,dimw-1 do begin
              dslice_cgw[*,*,j,i] = frebin(reform(slice_w[*,*,j,i]),nx_cg,ny_cg,/total)
           endfor
           print,string(13b)+' % finished: ',float(i)*100./((num_t-1)>1),format='(a,f4.0,$)'
        endfor
        for i=0,nx_cg-1 do begin
           for j=0,ny_cg-1 do begin
              dslice_cg[i,j,*,*] = frebin(reform(dslice_cgw[i,j,*,*]),nw_cg,num_t,/total)
           endfor
           print,string(13b)+' % finished: ',float(i)*100./((nx_cg-1)>1),format='(a,f4.0,$)'
        endfor
     endif else begin
        slice_s = slice_t*0.
        for i=0,num_t-1 do begin
           if max(slice_t[*,*,i]) gt 0 then begin
              slice_s[*,*,i] = filter_image(reform(slice_t[*,*,i]),FWHM=[numpix_x_fwhm,numpix_w],/all)
           endif
           print,string(13b)+' % finished: ',float(i)*100./((num_t-1)>1),format='(a,f4.0,$)'
        endfor
        for i=0,num_t-1 do dslice_cg[*,*,i] = frebin(reform(slice_s[*,*,i]),nx_cg,nw_cg,/total)
     endelse
     dslice_cg = dslice_cg*dw*effarea*apix
  endif else begin
     if dimy ne 0 then begin
        slice_s = slice_t*0.
        for i=0,num_t-1 do begin
           if max(slice_t[*,*,i]) gt 0 then begin
              slice_s[*,*,i] = filter_image(reform(slice_t[*,*,i]),FWHM=[numpix_x_fwhm,numpix_y_fwhm],/all)
           endif
           print,string(13b)+' % finished: ',float(i)*100./((num_t-1)>1),format='(a,f4.0,$)'
        endfor
        for i=0,num_t-1 do dslice_cg[*,*,i] = frebin(reform(slice_s[*,*,i]),nx_cg,ny_cg,/total)
     endif else begin
        slice_s = filter_image(slice_t,FWHM=[numpix_x_fwhm,1],/all)
        dslice_cg = frebin(slice_s,nx_cg,num_t,/total)
     endelse
     dslice_cg = dslice_cg*dw*effarea*apix
  endelse        

  if dimw ne 0 then int_dslice = total(dslice_cg,locw)/cf[nwave/2] else int_dslice = dslice_cg/cf
  if keyword_set(noise) then begin
     int_dslice_n = poidev(int_dslice)
     
     dslice_cgn = dslice_cg*0.
; add random fluctuations on the order of wfluc % to each wavelength position: 
     if dimw ne 0 then begin
        if dimy ne 0 then begin           
           for i=0,num_t-1 do begin
              for j=0,nx_cg-1 do begin
                 for k=0,ny_cg do begin
                    int = int_dslice_n[j,k,i]
                    if int gt 0. then begin
                       rnd_sp = randomn(777.+i*j+k,nw_cg)*wfluc*dslice_cg[j,k,*,i]+dslice_cg[j,k,*,i]
                       rnd_tot = total(rnd_sp)
                       dslice_cgn[j,k,*,i] = rnd_sp*int/rnd_tot
                    endif
                 endfor
              endfor
              print,string(13b)+' % finished: ',float(i)*100./((num_t-1)>1),format='(a,f4.0,$)'
           endfor
        endif else begin
           for i=0,num_t-1 do begin
              for j=0,nx_cg-1 do begin
                 int = int_dslice_n[j,i]
                 if int gt 0. then begin
                    rnd_sp = randomn(777.+i*j,nw_cg)*wfluc*dslice_cg[j,*,i]+dslice_cg[j,*,i]
                    rnd_tot = total(rnd_sp)
                    dslice_cgn[j,*,i] = rnd_sp*int/rnd_tot
                 endif
              endfor
              print,string(13b)+' % finished: ',float(i)*100./((num_t-1)>1),format='(a,f4.0,$)'
           endfor
        endelse
     endif else begin
        dslice_cgn = int_dslice_n
     endelse
  endif else begin
     if dimw eq 0 then dslice_cgn = int_dslice else dslice_cgn = dslice_cg/cf
  endelse

end
