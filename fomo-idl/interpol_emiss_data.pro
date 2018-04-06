
pro interpol_emiss_data,ne_s=ne_s,te=te,ion=ion,gotdir=gotdir,w0=w0,emission_goft=emission_goft,g_logte=g_logte,g_logne=g_logne,tstep=tstep,sav=sav,sdir=sdir,filenm=filenm,file_abund=file_abund,channel=channel,silent=silent

if keyword_set(ion) eq 0 then begin
   print,'interpol_emiss_data,ne_s=ne_s,te=te,ion=ion,gotdir=gotdir,w0=w0,emission_goft=emission_goft,g_logte=g_logte,g_logne=g_logne,tstep=tstep,sav=sav,sdir=sdir,filenm=filenm,file_abund=file_abund,channel=channelsilent=silent'
   return
endif

; Performs bicubic interpolation to given temperature and density
; points based on a given G(T,n) contribution function for a specific
; line emission (or imaging channel, e.g. AIA or DKIST filter), and
; returns emissivity function 

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
      nab = '_abph'
   endif
   if file_abund eq 'coronal' then begin
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
   nab = '_abco'
endelse

emission_goft = ne_s*0.d & interp_goft = ne_s*0.d
siz = size(ne_s)
dimx = siz[1] & dimy = siz[2]

lookup_goft,ion=ion,w0=w0,gotdir=gotdir,n_e_lg=g_logne,logt=g_logte,goft_mat=goft_mat,watom= watom,nab=nab,filenm=filenm,channel=channel,silent=silent

logne = alog10(ne_s)
logte = alog10(te)

;interpolation
goft_mat_size = size(goft_mat)
logne_ind = interpol(findgen(goft_mat_size(1)), g_logne, logne)
logte_ind = interpol(findgen(goft_mat_size(2)), g_logte, logte)
emission_goft =  interpolate(goft_mat,logne_ind,logte_ind, cubic = -0.5) * ne_s^2
;end of interpolation


if keyword_set(sav) and ~keyword_set(channel) then begin
   ntstep = string(tstep,"(i4)")
   if keyword_set(filenm) then fil1 = filenm else fil1 = ''
   if keyword_set(channel) then fil2 = filenm else fil2 = ''
   save,emission_goft,filename=sdir+'emission_goft'+fil1+'_'+fil2+'_'+ntstep+'.sav'
endif

end
