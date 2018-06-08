pro make_eitresponse, sngfilter=sngfilter, gotdir=gotdir, file_abund=file_abund,extname=extname

;if arg_present(sngfilter) lt 1 then begin
;   print,'make_eitresponse, sngfilter=sngfilter, gotdir=gotdir, file_abund=file_abund'
;   return
;endif

; Generates G(T,n) tables (200 x 200 pts) multiplied with EIT response functions and aperture size.
; Implies taking (dx/f)^2 in the further code (FoMo-C or FoMo-idl).
; It assumes the chianti.ioneq CHIANTI file for ionization equilibrium
; values. 

; INPUT:
; sngfilter = (string) 'all' or  name of filter (ex: '171' for the EIT 171 filter):
;              EUV:
;              '304' -> EIT 304
;              '171' -> EIT 171
;              '195' -> EIT 195
;              '284' -> EIT 284
; gotdir = (string) directory path where to save the generated table
;          (don't forget '/' at end of path) 
; file_abund = (string) 'coronal' or 'photospheric' depending on whether
;             'sun_coronal.abund' or 'sun_photospheric.abund' CHIANTI packages,
;             respectively, are to be used. 
; CALLS:
; eit_parms, isothermal

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

ioneq_name = concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq')

numt = 200
temp = 10.d^(findgen(numt)/(numt-1)*4.+4.0)
alogt = alog10(temp)

n_e_min = 1.e8
n_e_max = 1.e11

;steplg = 0.0015
steplg = 0.03
numn = alog10(n_e_max/n_e_min)/steplg

n_e_lg = dindgen(numn+1)/numn*alog10(n_e_max/n_e_min)+alog10(n_e_min)
n_e = 10.^(n_e_lg)

;sterad_eit_pix=(1.275e-5)^2 ;;;  EIT pixel - 2.629arcsec; sterad_eit_pix=(1.275e-5)^2=(2.629arcsec)^2
sterad_eit_pix=1.0
watom = 0.

NN=1000
wmin171=160.0 & wmax171=250.0 & lambda171=findgen(NN)/(NN-1)*(wmax171-wmin171)+wmin171
wmin195=160.0 & wmax195=250.0 & lambda195=findgen(NN)/(NN-1)*(wmax195-wmin195)+wmin195
wmin284=200.0 & wmax284=330.0 & lambda284=findgen(NN)/(NN-1)*(wmax284-wmin284)+wmin284
wmin304=200.0 & wmax304=350.0 & lambda304=findgen(NN)/(NN-1)*(wmax304-wmin304)+wmin304
eit_resp_171 = eit_parms(lambda171,171,units=units)
eit_resp_195 = eit_parms(lambda195,195,units=units)
eit_resp_284 = eit_parms(lambda284,284,units=units)
eit_resp_304 = eit_parms(lambda304,304,units=units)

if sngfilter eq 'all' or sngfilter eq '304' then begin
   openw,unit1,gotdir+'goft_table_eit304_'+nab+extname+'.dat',/get_lun & w0_1 = 304. & ion_1 = '304'
   printf,unit1,ion_1
   printf,unit1,w0_1
   printf,unit1,watom
   printf,unit1,numn,numt
   printf,unit1,alogt
endif
if sngfilter eq 'all' or sngfilter eq '171' then begin 
   openw,unit5,gotdir+'goft_table_eit171_'+nab+extname+'.dat',/get_lun & w0_5 = 171. & ion_5 = '171'
   printf,unit5,ion_5
   printf,unit5,w0_5
   printf,unit5,watom
   printf,unit5,numn,numt
   printf,unit5,alogt
endif
if sngfilter eq 'all' or sngfilter eq '195' then begin
   openw,unit6,gotdir+'goft_table_eit195_'+nab+extname+'.dat',/get_lun & w0_6 = 195. & ion_6 = '195'
   printf,unit6,ion_6
   printf,unit6,w0_6
   printf,unit6,watom
   printf,unit6,numn,numt
   printf,unit6,alogt
endif
if sngfilter eq 'all' or sngfilter eq '284' then begin
   openw,unit7,gotdir+'goft_table_eit284_'+nab+extname+'.dat',/get_lun & w0_7 = 284. & ion_7 = '284'
   printf,unit7,ion_7
   printf,unit7,w0_7
   printf,unit7,watom
   printf,unit7,numn,numt
   printf,unit7,alogt
endif

for i=0,numn-1 do begin

   ;isothermal, 160.0, 350.0, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file

   if sngfilter eq 'all' or sngfilter eq '304' then begin
      isothermal, wmin304, wmax304, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file
      eff_304=interpol(eit_resp_304,lambda304,lambda)
      ;eff_304=eit_resp_304
      sp_conv_304 = spectrum & sp_conv_304[*,*]=0. 
      for j=0,numt-1 do sp_conv_304[*,j]=fnhne*sterad_eit_pix*spectrum[*,j]*eff_304
      resp_304=total(sp_conv_304,1)
      printf,unit1,n_e_lg[i]
      printf,unit1,resp_304
   endif

   if sngfilter eq 'all' or sngfilter eq '171' then begin
      isothermal, wmin171, wmax171, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file
      eff_171=interpol(eit_resp_171,lambda171,lambda)
      ;eff_171=eit_resp_171
      sp_conv_171 = spectrum & sp_conv_171[*,*]=0. 
      for j=0,numt-1 do sp_conv_171[*,j]=fnhne*sterad_eit_pix*spectrum[*,j]*eff_171
      resp_171=total(sp_conv_171,1)
      printf,unit5,n_e_lg[i]
      printf,unit5,resp_171
   endif

   if sngfilter eq 'all' or sngfilter eq '195' then begin
      isothermal, wmin195, wmax195, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file
      eff_195=interpol(eit_resp_195,lambda195,lambda)
      ;eff_195=eit_resp_195
      sp_conv_195 = spectrum & sp_conv_195[*,*]=0. 
      for j=0,numt-1 do sp_conv_195[*,j]=fnhne*sterad_eit_pix*spectrum[*,j]*eff_195   
      resp_195=total(sp_conv_195,1)
      printf,unit6,n_e_lg[i]
      printf,unit6,resp_195
   endif    

   if sngfilter eq 'all' or sngfilter eq '284' then begin
      isothermal, wmin284, wmax284, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file
      eff_284=interpol(eit_resp_284,lambda284,lambda)
      ;eff_284=eit_resp_284
      sp_conv_284 = spectrum & sp_conv_284[*,*]=0.
      for j=0,numt-1 do sp_conv_284[*,j]=fnhne*sterad_eit_pix*spectrum[*,j]*eff_284
      resp_284=total(sp_conv_284,1)
      printf,unit7,n_e_lg[i]
      printf,unit7,resp_284
   endif
 
   print,string(13b)+' % finished: ',float(i)*100./(numn-1),format='(a,f4.0,$)'
endfor

if sngfilter eq 'all' or sngfilter eq '304' then free_lun,unit1
if sngfilter eq 'all' or sngfilter eq '171' then free_lun,unit5 
if sngfilter eq 'all' or sngfilter eq '195' then free_lun,unit6
if sngfilter eq 'all' or sngfilter eq '284' then free_lun,unit7

end
