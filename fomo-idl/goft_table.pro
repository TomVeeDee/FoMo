

pro goft_table,w0=w0,ion=ion,gotdir=gotdir,file_abund=file_abund,silent=silent,fact=fact,num=num

if keyword_set(w0) eq 0 or keyword_set(ion) eq 0 then begin
   print,'goft_table, w0=w0, ion=ion,gotdir=gotdir,file_abund=file_abund,silent=silent'
   return
endif

; Generates a G(T,n) table (200 pts in temperature, 3000 points in density) for a given line transition.
; A coronal abundance is included by default ('sun_coronal_2012_schmelz.abund').

; INPUT:
; w0 = (float) line center wavelength of line transition in Angstroms
; ion = (string) acronym of ion (see elements.pro)
; gotdir = (string) directory path where to save the generated table
;          (don't forget '/' at end of path) 
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
      nab = 'abph'
   endif
   if file_abund eq 'coronal' then begin
      abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_2012_schmelz.abund') ;!xuvtop+'/abundance/sun_coronal.abund'
      if file_test(abund_name) eq 0 then begin
         abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal_2012_schmelz.abund")'
      endelse
      nab = 'abco'
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
   abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_2012_schmelz.abund')
   if file_test(abund_name) eq 0 then begin
         abund_name = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal_2012_schmelz.abund")'
      endelse
      nab = 'abco'
endelse

; Ionisation file:
ioneq_name = concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq')

; determine wavelength transition:

;method1:
;elements,w0=w0,cw0=cw0,ion=ion,enum=enum,inum=inum,ind=ind,lv1=lv1,lv2=lv2

; method 2:

which_line_fomo,ion=ion,w0=w0,cw0=cw0,lv1=lv1,lv2=lv2,all=all,silent=silent,fact=fact,num=num,enum=enum

if w0 ne cw0 then begin
   print,'Warning: input wavelength (Angs.):',w0
   print,'Corresponding wavelength in Chianti (Angstrom):',cw0
endif else begin
   print,'Corresponding wavelength in Chianti (Angstrom):',cw0
endelse

ne_lg = 9
gofnt,ion,cw0-0.0001,cw0+0.0001,temp,goft0,desc,dens=10.^ne_lg,ioneq_name=ioneq_name,abund_name=abund_name,lower_levels=lv1,upper_levels=lv2,verbose=silent

;alogt0 = findgen(101)/20+4. ; same range as that defined by Chianti
alogt0 = alog10(temp)
lclgtm = ([max(goft0),!c])[1]
logTm = alogt0[lclgtm]

  ; get atomic weight:
watom = get_atomic_weight(enum)

n_e_min = 1.e8
if logtm lt 5 then n_e_max = 1.e12 else n_e_max = 1.e11

steplg = 0.001
;steplg = 0.002

units = 'erg cm^3 s^{-1}'
openr,unitversion, !xuvtop+'/VERSION',/get_lun
str=''
readf,unitversion, str
free_lun,unitversion
vchianti = 'CHIANTI'+str

numn = round(alog10(n_e_max/n_e_min)/steplg)

n_e_lg = dindgen(numn+1)/numn*alog10(n_e_max/n_e_min)+alog10(n_e_min)

if round(w0) lt 9999 then w0nm = string(round(w0),format='(i4.4)') else w0nm = string(round(w0),format='(i5.5)') 
openw,unit,gotdir+'goft_table_'+ion+'_'+w0nm+'_'+nab+'.dat',/get_lun
printf,unit,ion
printf,unit,cw0
printf,unit,watom
printf,unit,lv1,lv2
printf,unit,units
printf,unit,vchianti

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
   if ~keyword_set(silent) then print,'doing density ', i, ' of ',numn
   gofnt,ion,cw0-0.0001,cw0+0.0001,temp,goft0,desc,dens=10.^n_e_lg[i],ioneq_name=ioneq_name,abund_name=abund_name,lower_levels=lv1,upper_levels=lv2,verbose=silent

;   goft1 = goft0[pts]
;   ion_interp,alogt1,goft1,alogt2,goft2
   alogt = alog10(temp)
   goft2 = (interpol(goft0,alogt,alogt2,/spline))>0.

   printf,unit,n_e_lg[i]
   printf,unit,goft2
endfor
free_lun,unit

end
