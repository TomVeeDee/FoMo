pro make_dkistresponse, sngfilter=sngfilter, wvlmin=wvlmin, wvlmax=wvlmax, gotdir=gotdir, file_abund=file_abund,extname=extname,extro=extro,silent=silent

if ~keyword_set(sngfilter) then begin
   print,'make_dkistresponse, sngfilter=sngfilter, wvlmin=wvlmin, wvlmax=wvlmax, gotdir=gotdir, file_abund=file_abund'
   return
endif

; Generates G(T,n) tables (200 x 200 pts) for the DKIST response functions. 
; It assumes the chianti.ioneq CHIANTI file for ionization equilibrium
; values. 

; INPUT:
; sngfilter = (string) 'all' if all filters are to be generated

; '10747' -> Fe XIII 10747 A line (Filter FWHM: 10 A)
; '39340' -> Si IX 39340 A line (Filter FWHM: 200 A)

; wvlmin = (float) Minimum of desired wavelength range for line
; transition in Angstroms. Default: w0-5 for 10747 and w0-100 for 39340
; wvlmax = (float) Maximum of desired wavelength range for line
; transition in Angstroms. Default: w0+5 for 10747 and w0+100 for 39340
; gotdir = (string) directory path where to save the generated table
;          (don't forget '/' at end of path) 
; file_abund = (string) 'coronal' or 'photospheric' depending on whether
;             'sun_coronal.abund' or 'sun_photospheric.abund' CHIANTI packages,
;             respectively, are to be used. 
; extro: (string) for G(T,n) tables with extended density range [6,12]
;        in log. Default is [8,11] 

; CALLS:
; aia_get_response, isothermal

if ~keyword_set(extname) then extname = ''
if ~keyword_set(file_abund) then begin
   abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_2012_schmelz.abund')
   fnhne = 0.887 ; ratio of hydrogen to electrons number (= proton_dens(6.0) from Chianti, with coronal abundances)
   if file_test(abund_file) eq 0 then begin
      abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
      if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal.abund")'
   endif else begin
      if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal_2012_schmelz.abund")'
   endelse
   nab = 'abco'
endif else begin
   if file_abund eq 'coronal' then begin
      abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_2012_schmelz.abund') ;!xuvtop+'/abundance/sun_coronal.abund'
      fnhne = 0.887 ; ratio of protons to electrons number (= proton_dens(6.0) from Chianti, with coronal abundances)
      if file_test(abund_file) eq 0 then begin
         abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming coronal abundances (file:"sun_coronal_2012_schmelz.abund")'
      endelse
      nab = 'abco'
   endif
   if file_abund eq 'photospheric' then begin
      abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_photospheric.abund');!xuvtop+'/abundance/sun_photospheric.abund'
      fnhne = 0.848 ; ratio of protons to electrons number (= proton_dens(6.0) from Chianti, with photospheric abundances)
      if file_test(abund_file) eq 0 then begin
         abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_photospheric_1998_grevesse.abund')
         if ~keyword_set(silent) then print,'Assuming photospheric abundances (file:"sun_photospheric_1998_grevesse.abund")'
      endif else begin
         if ~keyword_set(silent) then print,'Assuming photospheric abundances (file:"sun_photospheric.abund")'
      endelse
      nab = 'abph'
   endif
endelse

; this is needed for the correct execution of CHIANTI's idl/utilities/proton_dens.pro
!abund_file = abund_file

ioneq_name = concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq')

numt = 200
temp = 10.d^(findgen(numt)/(numt-1)*4.+4.0)
alogt = alog10(temp)

if keyword_set(extro) then begin
   n_e_min = 1.e6
   n_e_max_sml = 1.e12
   extname = extname+'_extro'
endif else begin
   n_e_min = 1.e8
   n_e_max_sml = 1.e11
endelse

steplg = 0.015
numn_sml = round(alog10(n_e_max_sml/n_e_min)/steplg)

n_e_lg_sml = dindgen(numn_sml+1)/numn_sml*alog10(n_e_max_sml/n_e_min)+alog10(n_e_min)
n_e_sml = 10.^(n_e_lg_sml)   

; solid angle of 1 arcsec^2:
;sterad_arc=(!pi/180./3600.)^2

; normalised over the sphere:
sterad_arc = 1.
watom = 0.
units = 'cm^5 DN s^{-1} sr^{-1}'
openr,unitversion, !xuvtop+'/VERSION',/get_lun
str=''
readf,unitversion, str
free_lun,unitversion
vchianti = 'CHIANTI'+str

if sngfilter eq 'all' then dkarr = ['10747','39340'] else dkarr = sngfilter
ndkar = n_elements(dkarr)

for k=0,ndkar-1 do begin
   filt = dkarr[k]
   print,'Constructing DKIST - '+filt+' G(T) table'
   openw,unit,gotdir+'goft_table_dkist'+filt+'_'+nab+extname+'.dat',/get_lun & w0 = float(filt) & ion = filt
   if ~keyword_set(wvlmin) or ~keyword_set(wvlmax) or keyword_set(all) then begin
      if filt eq '10747' then begin
         wvlmin = w0-5.
         wvlmax = w0+5.
      endif
      if filt eq '39340' then begin
         wvlmin = w0-100.
         wvlmax = w0+100.
      endif
   endif
   printf,unit,ion
   printf,unit,w0
   printf,unit,watom
   printf,unit,wvlmin,wvlmax
   printf,unit,units
   printf,unit,vchianti
   printf,unit,numn_sml,numt
   printf,unit,alogt

   for i=0,numn_sml-1 do begin
      print,'DKIST '+filt+' - doing density ',i,' of ',numn_sml
      isothermal, wvlmin, wvlmax, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e_sml[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file,noverbose=silent
      sp_conv = spectrum & sp_conv[*,*]=0.
      for j=0,numt-1 do sp_conv[*,j] = fnhne*sterad_arc*spectrum[*,j]
      resp = total(sp_conv,1)
      printf,unit,n_e_lg_sml[i]
      printf,unit,resp
   endfor
   free_lun,unit
endfor

end
