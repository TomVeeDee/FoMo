

pro goft_table,w0=w0,ion=ion,gotdir=gotdir,file_abund=file_abund,vers=vers,silent=silent

if keyword_set(w0) eq 0 or keyword_set(ion) eq 0 then begin
   print,'goft_table, w0=w0, ion=ion,gotdir=gotdir,file_abund=file_abund,vers=vers,silent=silent'
   return
endif

; Generates a G(T,n) table (200 pts in temperature, 3000 points in density) for a given line transition.
; A coronal abundance is included by default ('sun_coronal.abund').

; INPUT:
; w0 = (float) line center wavelength of line transition in Angstroms
; ion = (string) acronym of ion (see elements.pro)
; gotdir = (string) directory path where to save the generated table
;          (don't forget '/' at end of path) 
; vers = (single integer 6, 7 or 8). Chianti version. Read from
; installed Chianti package if not known
; file_abund: (string) file for abundance abundance. 2 kinds are implemented:
;            'photospheric' or 'coronal' corresponding, respectively, to the
;            CHIANTI packages: sun_coronal.abund and sun_photospheric.abund
;            By default the 'coronal' abundance package is set.
;	     If other abundance is desired set file_abund = 'other' and provide 
;	     the full path to the abundance file in the keyword 'ext_abund'.

; CALLS:
; elements, get_atomic_weight, g_of_t

if keyword_set(file_abund) then begin
   if file_abund eq 'photospheric' then begin
      abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_photospheric.abund');!xuvtop+'/abundance/sun_photospheric.abund'
      if file_test(abund_name) eq 0 then begin
         abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_photospheric_1998_grevesse.abund')
         if ~keyword_set(silent) then print,'Assuming photospheric abundances (file:"sun_photospheric_1998_grevesse.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming photospheric abundances (file:"sun_photospheric.abund")'
      endelse
      nab = '_abph'
   endif
   if file_abund eq 'coronal' then begin
      abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund') ;!xuvtop+'/abundance/sun_coronal.abund'
      if file_test(abund_name) eq 0 then begin
         abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_2012_schmelz.abund')
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal_1992_feldman.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal.abund")'
      endelse
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
   abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
   if file_test(abund_name) eq 0 then begin
         abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_2012_schmelz.abund')
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal_1992_feldman.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal.abund")'
      endelse
      nab = 'abco'
endelse

  ; determine wavelength transition:
elements, w0=w0, ion=ion, logTm=logTm, enum=enum, inum=inum, ind=ind, vers=vers

if ~keyword_set(ind) then begin
   emiss = emiss_calc(enum,inum)
   n_elt = n_elements(emiss.lambda)
   index = lindgen(n_elt)
   wavels = emiss(index).lambda
   ind = ([min(abs(wavels-w0)),!c])[1]
endif
  ; get atomic weight:
  watom = get_atomic_weight(enum)

  n_e_min = 1.e8
;  n_e_max = 1.e11
  n_e_max = 1.e12

;steplg = 0.001
steplg = 0.002

numn = alog10(n_e_max/n_e_min)/steplg

n_e_lg = dindgen(numn+1)/numn*alog10(n_e_max/n_e_min)+alog10(n_e_min)

w0nm = string(round(w0),format='(i4.4)')
openw,unit,gotdir+'goft_table_'+ion+'_'+w0nm+'_'+nab+'.dat',/get_lun
printf,unit,ion
printf,unit,w0
printf,unit,watom

;tmin = (logTm-wte)>4.           ;define wanted range
wte = 0.75
if logTm-wte lt 4 then tmin = 3.5 else tmin = (logTm-wte)>4.
alogt0 = findgen(101)/20+4. ; same range as that defined by Chianti
tmax = (logTm+wte)<7.
pts = where(alogt0 ge tmin and alogt0 le tmax,npts)
;alogt1 = alogt0[pts]
numt = 200
alogt2 = findgen(numt)/(numt-1)*2*wte+tmin

printf,unit,numn,numt
printf,unit,alogt2
for i=0,numn do begin
   goft0=g_of_t(enum,inum,dens=n_e_lg[i],ioneq_file=concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq'),abund_file=abund_name,index=ind,/quiet)

;   goft1 = goft0[pts]
;   ion_interp,alogt1,goft1,alogt2,goft2
   goft2 = (interpol(goft0,alogt0,alogt2))>0.

   printf,unit,n_e_lg[i]
   printf,unit,goft2
endfor
free_lun,unit

end
