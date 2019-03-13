pro read_line, dir=dir,idcase=idcase,w0=w0,ion=ion,gridx=gridx,gridy=gridy,line1_t=line1_t,line2_t=line2_t,line3_t=line3_t,time=time,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,grid_1=grid_1,grid_2=grid_2,grid_3=grid_3,npixlos_1=npixlos_1,npixlos_2=npixlos_2,npixlos_3=npixlos_3,wave=wave,mua_d=mua_d,emission_goft_t=emission_goft_t,channel=channel,imaging=imaging,lcx=lcx,lct=lct,box=box

if ~keyword_set(dir) then begin
   print,'read_line, dir=dir,idcase=idcase,w0=w0,ion=ion,gridx=gridx,gridy=gridy,line1_t=line1_t,line2_t=line2_t,line3_t=line3_t,time=time,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,grid_1=grid_1,grid_2=grid_2,grid_3=grid_3,npixlos_1=npixlos_1,npixlos_2=npixlos_2,npixlos_3=npixlos_3,wave=wave,mua_d=mua_d,emission_goft_t=emission_goft_t,channel=channel,imaging=imaging,lcx=lcx,lct=lct,box=box'
   return
endif

; INPUT
; this routine reads from files params_* and fomo_* files that are
; generated from fomo_synth calculation. 

; dir:  input data directory

; idcase: name of specific case being analysed (assumes an existing
; directory with this name under directory 'dir')

; w0: rest wavelength of spectral line transition

; ion: ion of spectral line 

; OUTPUT:

; gridx: x-grid of original numerical box

; gridy: y-grid of original numerical box

; line1_t: the specific intensity, result of fomo_synth
; calculation along the LOS mua_d[0]. 
; If the calculation invoves spectra then this is a 
; (dim_1,nwave,nstep) array, where: 
; dim_1 is the dimension perpendicular to the line-of-sight (LOS)
; nwave is the number of elements in wavelength
; nstep is the number of time steps
; If the calculation is only imaging then line1_t is (dim_1,nstep)

; Similarly for line2_t and line3_t where the index number denotes
; different LOS angles mua_d[1] and mua_d[2] (if any). 
; If more than 3 angles in the calculation,
; these need to be explicitly added as keywords to the output.

; mua_d: the array with LOS angles considered in the fomo_synth run.

; time: time array in s

; dl_1: the spatial unit (CGS) in the uniform grid_1 (same for dl_2
; and dl_3 for grid_2 and grid_3)

; grid_1: the grid perpendicular to the LOS mua_d[0] 
; (same for grid_2 and grid_3 for the corresponding LOS angles)

; wave: the wavelength array if the calculation has spectra (no imaging)

; emission_goft_t: the emissivity array G(T)*ne^2 used in the
; calculation of the line intensities. This has the same units 
; as the numerical box (nx,ny,nstep)

; set 'channel' keyword to the name of the channel used when
; performing imaging calculations (no spectra)

; set 'imaging' keyword when reading from imaging calculations (not
; necessary if 'channel' is set. 'imaging' is needed when no specific
; broadband channel is modelled (like AIA 171) but only a specific
; line transition without spectra (so total intensity only).

; OPTIONAL:
; /box: due to the numerical box limitation, the integration along
; oblique LOS has different integration lengths, which results in
; different emission from the edge to the centre of the box. By
; setting this keyword this effect is corrected for by making all
; paths the same length. This is made by assuming that the pixel(s)
; crossed by a LOS ray at the very edge of the domain correspond to
; emission at rest. The emissivity at rest from a single pixel is 
; calculated and added to each LOS ray a number of time equal to match
; the longest LOS ray. The number of times for each LOS ray is
; returned in 'pixlos_1' (respectively, pixlos_2 and pixlos_3 for the
; other LOS angles). 

; lct: a specific time step to read from
; lcx: a specific z-cut to read from 

savdir = dir+idcase+'/sav/'
if (keyword_set(imaging) or keyword_set(channel)) then name = '_'+idcase+'_imag' else name = '_'+idcase+'_synth'

if w0 lt 1.e3 then w0nm = string(w0,format="(i3.3)") else w0nm = string(w0,format="(i4.4)")
if w0 lt 1.e4 and ~keyword_set(channel) then w0nm = string(round(w0),format='(i4.4)') else w0nm = string(round(w0),format='(i3.3)')
if w0 gt 1.e4 then w0nm = string(round(w0),format='(i5.5)')
nw0 = string(w0,format='(d0.3)')

if keyword_set(channel) or keyword_set(imaging) then begin
   if keyword_set(channel) then name = name+'_'+channel+'_'+w0nm else name = name +'_'+ion+'_'+nw0
endif else begin
   name = name +'_'+ion+'_'+nw0
endelse

restore,savdir+'params_fomo'+name+'.sav'

if n_elements(lcx) then nlcx = '_'+string(lcx,format="(i5.5)")
if n_elements(lct) then begin
   nlct = string(lct,format="(i4.4)")
   print,'Reading values for time = '+nlct
endif
if n_elements(lct) then files = file_search(savdir+'fomo'+name+'*_'+nlct+'.sav',count=nfiles,/fully_qualify_path) 

if n_elements(lcx) then files = file_search(savdir+'fomo'+name+'*_'+nlcx+'.sav',count=nfiles,/fully_qualify_path)

if n_elements(lcx) eq 0 and n_elements(lct) eq 0 then files = file_search(savdir+'fomo'+name+'*.sav',count=nfiles,/fully_qualify_path)

restore,files[-1]

nmua = n_elements(mua_d)

if ~keyword_set(imaging) and ~keyword_set(channel) then nwave = (size(line_1))[2] else nwave = 1

if n_elements(lct) then nstep = 1 else nstep = n_elements(files)

time = fltarr(nstep)

if ~keyword_set(imaging) and ~keyword_set(channel) then begin
   for j=0,nmua-1 do begin
      sn = string(j+1,format="(i1)") 
      ex1 = 'dim_'+sn+' = (size(line_'+sn+'))[1]'
      void = execute(ex1)
      ex1 = 'line'+sn+'_t = fltarr(dim_'+sn+',nwave,nstep)'
      void = execute(ex1)
   endfor
endif else begin
   for j=0,nmua-1 do begin
      sn = string(j+1,format="(i1)") 
      ex1 = 'dim_'+sn+' = (size(line_'+sn+'))[1]'
      void = execute(ex1)
      ex1 = 'line'+sn+'_t = fltarr(dim_'+sn+',nstep)'
      void = execute(ex1)
   endfor
endelse
emission_goft_t = fltarr(n_elements(gridx),n_elements(gridy),nstep)

for i=0,nstep-1 do begin
   restore,files[i]
   time[i] = tstep
   if ~keyword_set(imaging) and ~keyword_set(channel) then begin
      for j=0,nmua-1 do begin
         sn = string(j+1,format="(i1)") 
         ex1 = 'npixlos_'+sn+' = intarr(dim_'+sn+')'
         void = execute(ex1)
         ex1 = 'line'+sn+'_t[*,*,i] = line_'+sn
         void = execute(ex1)
         ex1 = 'grid_'+sn+' = dl_'+sn+'*findgen(dim_'+sn+')'
         void = execute(ex1)
         if keyword_set(box) then begin
            ex1 = 'npixlos = npixlos_'+sn
            void = execute(ex1)
            ex1 = 'line = line_'+sn
            void = execute(ex1)
            ex1 = 'dim_ = dim_'+sn
            void = execute(ex1)
            ex1 = 'for k=0,dim_'+sn+'-2 do npixlos[k] = n_elements(emission_goft[n_gridx_'+sn+'[ngrid_'+sn+'[k]:ngrid_'+sn+'[k+1]-1],n_gridy_'+sn+'[ngrid_'+sn+'[k]:ngrid_'+sn+'[k+1]-1]])'
            void = execute(ex1)
            dm = reform(line[0,*]/npixlos[0])
            mxnpix = max(npixlos)
            if mxnpix ne npixlos[0] then begin
               for k=0,dim_-1 do begin
                  hm_npix_ray = mxnpix-npixlos[k]
                  ex1 = 'line'+sn+'_t[k,*,i] = line'+sn+'_t[k,*,i]+dm*total(hm_npix_ray)'
                  void = execute(ex1)
               endfor
            endif
         endif
      endfor
   endif else begin
      for j=0,nmua-1 do begin
         sn = string(j+1,format="(i1)") 
         ex1 = 'npixlos_'+sn+' = intarr(dim_'+sn+')'
         void = execute(ex1)
         ex1 = 'line'+sn+'_t[*,i] = line_'+sn
         void = execute(ex1)
         ex1 = 'grid_'+sn+' = dl_'+sn+'*findgen(dim_'+sn+')'
         void = execute(ex1)
         if keyword_set(box) then begin
            ex1 = 'npixlos = npixlos_'+sn
            void = execute(ex1)
            ex1 = 'line = line_'+sn
            void = execute(ex1)
            ex1 = 'dim_ = dim_'+sn
            void = execute(ex1)
            ex1 = 'for k=0,dim_'+sn+'-2 do npixlos[k] = n_elements(emission_goft[n_gridx_'+sn+'[ngrid_'+sn+'[k]:ngrid_'+sn+'[k+1]-1],n_gridy_'+sn+'[ngrid_'+sn+'[k]:ngrid_'+sn+'[k+1]-1]])'
            void = execute(ex1)
            dm = line[0]/npixlos[0]
            mxnpix = max(npixlos)
            if mxnpix ne npixlos[0] then begin
               for k=0,dim_-1 do begin
                  hm_npix_ray = mxnpix-npixlos[k]
                  ex1 = 'line'+sn+'_t[k,i] = line'+sn+'_t[k,i]+dm*total(hm_npix_ray)'
                  void = execute(ex1)
               endfor
            endif
         endif
      endfor
   endelse
   emission_goft_t[*,*,i] = emission_goft
   print,string(13b)+' % finished: ',float(i)*100./((nstep-1)>1),format='(a,f4.0,$)'
endfor 

end

