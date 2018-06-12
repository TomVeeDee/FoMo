pro make_aiaresponse, sngfilter=sngfilter, wvlmin=wvlmin, wvlmax=wvlmax, gotdir=gotdir, file_abund=file_abund,extname=extname

if ~keyword_set(sngfilter) then begin
   print,'make_aiaresponse, sngfilter=sngfilter, wvlmin=wvlmin, wvlmax=wvlmax, gotdir=gotdir, file_abund=file_abund'
   return
endif

; Generates G(T,n) tables for the AIA response functions (200 x 201
; for coronal channels or 200x267 for uv or 304 channels). 
; It assumes the chianti.ioneq CHIANTI file for ionization equilibrium
; values. 
; units: cm^5 DN s^-1 sr^-1 (no calculation by solid angle of 1 arcsec^2)

; INPUT:
; sngfilter = (string) 'all' if all filters, including (4500, 1700, 1600)
;              name of filter (ex: '171' for the AIA 171 filter):
;              EUV:
;              '304' -> AIA 304
;              '171' -> AIA 171
;              '193' -> AIA 193
;              '211' -> AIA 211
;              '335' -> AIA 335
;              '094' -> AIA 094
;              '131' -> AIA 131
;              (UV:)
;              '1600' -> AIA 1600
;              '1700' -> AIA 1700
;              '4500' -> AIA 4500
; wvlmin = (float) Minimum of desired wavelength range for single channel
; in Angstroms ; default: Wavelength value at 1% of peak
; wvlmax = (float) Maximum of desired wavelength range for single channel
; in Angstroms; default: Wavelength value at 1% of peak
;    note: wvlmin and wvlmax are calculated as above for sngfilter='all'
; gotdir = (string) directory path where to save the generated table
;          (don't forget '/' at end of path) 
; file_abund = (string) 'coronal' or 'photospheric' depending on whether
;             'sun_coronal.abund' or 'sun_photospheric.abund' CHIANTI packages,
;             respectively, are to be used. 
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

n_e_min = 1.e8
n_e_max_sml = 1.e11
n_e_max_big = 1.e12

steplg = 0.015
numn_sml = alog10(n_e_max_sml/n_e_min)/steplg
numn_big = alog10(n_e_max_big/n_e_min)/steplg

n_e_lg_sml = dindgen(numn_sml+1)/numn_sml*alog10(n_e_max_sml/n_e_min)+alog10(n_e_min)
n_e_lg_big = dindgen(numn_big+1)/numn_big*alog10(n_e_max_big/n_e_min)+alog10(n_e_min)
n_e_sml = 10.^(n_e_lg_sml)   
n_e_big = 10.^(n_e_lg_big)   

; solid angle of 1 arcsec^2:
;sterad_arc=(!pi/180./3600.)^2

; normalised over the sphere:
sterad_arc = 1.
watom = 0.

if sngfilter eq 'all' or sngfilter eq '1600' or sngfilter eq '1700' or sngfilter eq '4500' then aia_resp = aia_get_response(/dn,/uv) else aia_resp = aia_get_response(/dn)

if sngfilter eq 'all' then aiarr = ['304','1600','1700','4500','171','193','211','335','094','131'] else aiarr = sngfilter
naiar = n_elements(aiarr)

for k=0,naiar-1 do begin
   filt = aiarr[k]
   print,'Constructing AIA - '+filt+' G(T) table'

   openw,unit,gotdir+'goft_table_aia'+filt+'_'+nab+extname+'.dat',/get_lun & w0 = float(filt) & ion = filt
   if ~keyword_set(wvlmin) or ~keyword_set(wvlmax) or keyword_set(all) then begin
      if filt eq '094' then filt2 = '94' else filt2 = filt
      ex1 = 'pk = max(aia_resp.a'+filt2+'.ea,ipk)'
      ex2 = 'wave = aia_resp.a'+filt2+'.wave'
      void = execute(ex1)
      void = execute(ex2)
      ex1 = 'wvlmin = wave[([min(abs(pk*0.01-aia_resp.a'+filt2+'.ea[0:ipk])),!c])[1]]'
      ex2 = 'wvlmax = wave[([min(abs(pk*0.01-aia_resp.a'+filt2+'.ea[ipk:-1])),!c])[1]+ipk]'
   endif
   if filt eq '304' or filt eq '1600' or filt eq '1700' or filt eq '4500' then begin 
      numn = numn_big 
      n_e = n_e_big
      n_e_lg = n_e_lg_big
   endif else begin
      numn = numn_sml
      n_e = n_e_sml
      n_e_lg = n_e_lg_sml
   endelse
   printf,unit,ion
   printf,unit,w0
   printf,unit,watom
   printf,unit,wvlmin
   printf,unit,wvlmax
   printf,unit,numn,numt
   printf,unit,alogt

   for i=0,numn-1 do begin
      print,'AIA '+filt+' - doing density number ',i,' of ',numn
      isothermal, wvlmin, wvlmax, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file
      ex = 'eff = interpol(aia_resp.a'+filt2+'.ea,aia_resp.a'+filt2+'.wave,lambda,/spline)'
      void = execute(ex)
      sp_conv = spectrum & sp_conv[*,*]=0.
      for j=0,numt-1 do sp_conv[*,j] = fnhne*sterad_arc*spectrum[*,j]*eff
      resp = total(sp_conv,1)
      printf,unit,n_e_lg[i]
      printf,unit,resp
   endfor
   free_lun,unit
endfor

end
