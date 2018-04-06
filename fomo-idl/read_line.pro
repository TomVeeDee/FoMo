pro read_line, dir=dir,idcase=idcase,w0=w0,ion=ion,gridx=gridx,gridy=gridy,line1_t=line1_t,line2_t=line2_t,line3_t=line3_t,time=time,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,grid_1=grid_1,grid_2=grid_2,grid_3=grid_3,wave=wave,mua_d=mua_d,emission_goft_t=emission_goft_t,channel=channel,imaging=imaging,lcx3=lcx3,lct=lct

if ~keyword_set(dir) then begin
   print,'read_line, dir=dir,idcase=idcase,w0=w0,ion=ion,gridx=gridx,gridz=gridz,line1_t=line1_t,line2_t=line2_t,line3_t=line3_t,time=time,dl_1=dl_1,dl_2=dl_2,dl_3=dl_3,grid_1=grid_1,grid_2=grid_2,grid_3=grid_3,wave=wave,mua_d=mua_d,emission_goft_t=emission_goft_t,channel=channel,imaging=imaging,lcx3=lcx3,lct=lct'
   return
endif

; dir = input data directory

savdir = dir+idcase+'/sav/'
if (keyword_set(imaging) or keyword_set(channel)) then name = '_'+idcase+'_imag' else name = '_'+idcase+'_synth'

if w0 lt 1.e3 then w0nm = string(w0,format="(i3.3)") else w0nm = string(w0,format="(i4.4)")
;if w0 lt 1.e4 and ~keyword_set(channel) then w0nm = string(round(w0),format='(i4.4)') else w0nm = string(round(w0),format='(i3.3)')
;if w0 gt 1.e4 then w0nm = string(round(w0),format='(i5)')
nw0 = string(w0,format='(d0.3)')

if keyword_set(channel) or keyword_set(imaging) then begin
   if keyword_set(channel) then name = name+'_'+channel+'_'+w0nm else name = name +'_'+ion+'_'+nw0
endif else begin
   name = name +'_'+ion+'_'+nw0
endelse

restore,savdir+'params_fomo'+name+'.sav'

if n_elements(lcx3) then nlcx3 = '_'+string(lcx3,format="(i3.3)") else nlcx3 = ''
if n_elements(lct) then begin
   nlct = string(lct,format="(i4.4)")
   print,'Reading values for time = '+nlct
endif
if n_elements(lct) then files = file_search(savdir+'fomo'+name+'*_'+nlct+'.sav',count=nfiles,/fully_qualify_path) else files = file_search(savdir+'fomo'+name+nlcx3+'*.sav',count=nfiles,/fully_qualify_path)

restore,files[-1]

dim_1 = (size(line_1))[1]
if n_elements(mua_d) ge 2 then dim_2 = (size(line_2))[1]
if n_elements(mua_d) ge 3 then dim_3 = (size(line_3))[1] 
if ~keyword_set(imaging) and ~keyword_set(channel) then nwave = (size(line_1))[2] else nwave = 1
nt=2
if n_elements(lct) then nstep = n_elements(files) else nstep = nt

time = fltarr(nt)

if ~keyword_set(imaging) and ~keyword_set(channel) then begin
   line1_t = fltarr(dim_1,nwave,nstep)
   if n_elements(mua_d) ge 2 then line2_t = fltarr(dim_2,nwave,nstep)
   if n_elements(mua_d) ge 3 then line3_t = fltarr(dim_3,nwave,nstep)
endif else begin
   line1_t = fltarr(dim_1,nstep)
   if n_elements(mua_d) ge 2 then line2_t = fltarr(dim_2,nstep)
   if n_elements(mua_d) ge 3 then line3_t = fltarr(dim_3,nstep)
endelse
emission_goft_t = fltarr(dim_1,dim_3,nstep)

for i=0,nstep-1 do begin
   restore,files[i]
   time[i] = tstep
   if ~keyword_set(imaging) and ~keyword_set(channel) then begin
      line1_t[*,*,i] = line_1
      if n_elements(mua_d) ge 2 then line2_t[*,*,i] = line_2
      if n_elements(mua_d) ge 3 then line3_t[*,*,i] = line_3
   endif else begin
      line1_t[*,i] = line_1
      if n_elements(mua_d) ge 2 then line2_t[*,i] = line_2
      if n_elements(mua_d) ge 3 then line3_t[*,i] = line_3
   endelse
   emission_goft_t[*,*,i] = emission_goft
   print,string(13b)+' % finished: ',float(i)*100./((nstep-1)>1),format='(a,f4.0,$)'
endfor

grid_1 = dl_1*findgen(dim_1)
if n_elements(mua_d) ge 2 then grid_2 = dl_2*findgen(dim_2)
if n_elements(mua_d) ge 3 then grid_3 = dl_3*findgen(dim_3)

end

