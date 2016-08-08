
pro interpol_emiss_data,ne_s=ne_s,te=te,ion=ion,gotdir=gotdir,w0=w0,emission_goft=emission_goft,g_logte=g_logte,g_logne=g_logne,tstep=tstep,sav=sav,sdir=sdir,filenm=filenm,file_abund=file_abund,aia=aia,silent=silent

if keyword_set(ion) eq 0 then begin
   print,'interpol_emiss_data,ne_s=ne_s,te=te,ion=ion,gotdir=gotdir,w0=w0,emission_goft=emission_goft,g_logte=g_logte,g_logne=g_logne,tstep=tstep,sav=sav,sdir=sdir,filenm=filenm,file_abund=file_abund,aia=aia,silent=silent'
   return
endif

; Performs bilinear interpolation to given temperature and density
; points based on a given G(T,n) contribution function for a specific
; line transition (or AIA filter), and returns emissivity function (G(t,n)*ne^2). 

; INPUT: 
; ne_s, t_e: (0-2d float arrays) density and temperature values where
;           to interpolate (in CGS)
; ion: (string) acronym of the ion
; w0: (float) wavelength of line center 

; OUTPUT:
; emission_goft: (0-2d float arrays) where interpolated values are returned
; g_logte & g_logne: (1d float arrays) the log(T) and log(n) values
;           where the G(T,n) function is defined. 

; OPTIONAL:
; sav: set it for saving the calculated emissivity values.
; filenm & tstep: (string & float, resp.) used for naming the
;           emissivity sav file.  
; sdir: (string) path to where to save the emissivity table. 

; CALLS:
; lookup_goft

if keyword_set(file_abund) then begin
   if file_abund eq 'photospheric' then begin
      ; AIA goft tables calculated with abundance package: sun_photospheric_1998_grevesse.abund 
      nab = '_abph'
   endif
   if file_abund eq 'coronal' then begin
      ; AIA goft tables calculated with abundance package: sun_coronal_1992_feldman.abund
      nab = '_abco'
   endif
   if file_abund eq 'other' then begin
      if ~keyword_set(ext_abund) then begin
         print,'Please provide the name and full path of the abundance file in the keyword "ext_abund". Also, modify the extension "nab" accordingly'
         stop
      endif else begin
         abund_name = ext_abund
         nab = '_abext'         
      endelse
   endif
endif else begin
   ; Assuming coronal abundances
   ; AIA goft tables calculated with abundance package: sun_coronal_1992_feldman.abund
   nab = '_abco'
endelse

emission_goft = ne_s*0.d & interp_goft = ne_s*0.d
siz = size(ne_s)
dimx = siz[1] & dimy = siz[2]

lookup_goft,ion=ion,w0=w0,gotdir=gotdir,n_e_lg=g_logne,logt=g_logte,goft_mat=goft_mat,watom= watom,nab=nab,filenm=filenm,aia=aia,silent=silent

logne = alog10(ne_s)
logte = alog10(te)
lgne_sort = sort(logne)
logne_sorted = logne[lgne_sort] 
logte_sorted = logte[lgne_sort] 
arr_ne = logne_sorted*0. & num_arrne = n_elements(arr_ne)
arr_gne = findgen(n_elements(g_logne))
arr_te = logte_sorted*0. & num_arrte = n_elements(arr_te)
arr_gte = findgen(n_elements(g_logte))
for i=0,num_arrne-1 do arr_ne[lgne_sort[i]] = interpol(arr_gne,g_logne,logne_sorted[i])
for i=0,num_arrte-1 do arr_te[lgne_sort[i]] = interpol(arr_gte,g_logte,logte_sorted[i])

for i=0,num_arrne-1 do begin
   interp_goft[lgne_sort[i]] = interpolate(goft_mat,arr_ne[lgne_sort[i]],arr_te[lgne_sort[i]],/grid)
   print,string(13b)+' % finished: ',float(i)*100./(num_arrne-1),format='(a,f4.0,$)'
endfor
emission_goft[lgne_sort] = interp_goft[lgne_sort] * ne_s[lgne_sort]^2

if keyword_set(sav) eq 1 then begin
   ntstep = string(tstep,"(i4)")
   save,emission_goft,filename=sdir+'emission_goft'+filenm+'_'+ntstep+'.sav'
endif

end
