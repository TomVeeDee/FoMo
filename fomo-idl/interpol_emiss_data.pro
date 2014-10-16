
pro interpol_emiss_data,n_e,te,ion=ion, w0=w0,emission_goft=emission_goft, g_logte=g_logte, g_logne=g_logne, tstep=tstep, sav=sav, sdir=sdir, filenm=filenm,file_abund=file_abund,ext_abund=ext_abund

if keyword_set(ion) eq 0 then begin
   print,'interpol_emiss_data,n_e,te,ion=ion, w0=w0,emission_goft=emission_goft, g_logte=g_logte, g_logne=g_logne, tstep=tstep, sav=sav, sdir=sdir, filenm=filenm,file_abund=file_abund,ext_abund=ext_abund'
   return
endif

; Performs bilinear interpolation to given temperature and density
; points based on a given G(T,n) contribution function for a specific
; line transition (or AIA filter), and returns emissivity function (G(t,n)*ne^2)*abundance. 

; INPUT: 
; n_e, t_e: (0-2d float arrays) density and temperature values where
;           to interpolate (in CGS)
; ion: (string) acronym of the ion
; w0: (float) wavelength of line center 
; file_abund: (string) file for abundance abundance. 2 kinds are implemented:
;            'photospheric' or 'coronal' corresponding, respectively, to the
;            CHIANTI packages: sun_coronal.abund and sun_photospheric.abund
;            By default the 'coronal' abundance package is set.
;	     If other abundance is desired set file_abund = 'other' and provide 
;	     the full path to the abundance file in the keyword 'ext_abund'.

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
; lookup_goft, read_abund

if keyword_set(file_abund) then begin
   if file_abund eq 'photospheric' then begin
      abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_photospheric.abund');!xuvtop+'/abundance/sun_photospheric.abund'
      print,'Assuming photospheric abundances'
   endif
   if file_abund eq 'coronal' then begin
      abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund') ;!xuvtop+'/abundance/sun_coronal.abund'
      print,'Assuming coronal abundances'
   endif
   if file_abund eq 'other' then begin
      if ~keyword_set(ext_abund) then begin
         print,'Please provide the name and full path of the abundance file in the keyword "ext_abund"'
      endif else begin
         abund_name = ext_abund
      endelse
   endif
endif else begin
   abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
   print,'Assuming coronal abundances (file:"sun_coronal.abund")'
endelse

abund_dflt = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')

read_abund,abund_name,abund,abund_ref
line_abunds = abund[enum-1]
read_abund,abund_dflt,ab_dflt,abund_ref_dflt
line_abunds_dflt = ab_dflt[enum-1]

emission_goft = n_e*0.d & interp_goft = n_e*0.d
siz = size(n_e)
dimx = siz[1] & dimy = siz[2]

lookup_goft, ion=ion, w0=w0, n_e_lg=g_logne, logt=g_logte, goft_mat=goft_mat, watom= watom, filenm=filenm

logne = alog10(n_e)
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
emission_goft[lgne_sort] = interp_goft[lgne_sort] * n_e[lgne_sort]^2

if keyword_set(sav) eq 1 then begin
   ntstep = string(tstep,"(i4)")
   save,emission_goft,filename=sdir+'emission_goft'+filenm+'_'+ntstep+'.sav'
endif

emission_goft = emission_goft/line_abunds_dflt*line_abunds

end
