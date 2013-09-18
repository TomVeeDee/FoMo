
pro interpol_emiss_data,n_e,te,ion=ion, w0=w0,emission_goft=emission_goft, g_logte=g_logte, g_logne=g_logne, tstep=tstep, sav=sav, sdir=sdir, filenm=filenm

if keyword_set(ion) eq 0 then begin
   print,'interpol_emiss_data,n_e,te,ion=ion, w0=w0,emission_goft=emission_goft, g_logte=g_logte, g_logne=g_logne, tstep=tstep, sav=sav, sdir=sdir, filenm=filenm'
   return
endif

; Performs bilinear interpolation to given temperature and density
; points based on a given G(T,n) contribution function for a specific
; line transition (or AIA filter), and returns emissivity function (G(t,n)*ne^2). 

; INPUT: 
; n_e, t_e: (0-2d float arrays) density and temperature values where
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
; lookup_goft, 

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

end
