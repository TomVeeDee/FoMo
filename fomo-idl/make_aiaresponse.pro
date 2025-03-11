pro make_aiaresponse, sngfilter=sngfilter, wvlmin=wvlmin, wvlmax=wvlmax, gotdir=gotdir, file_abund=file_abund,extname=extname

if ~keyword_set(sngfilter) then begin
   print,'make_aiaresponse, sngfilter=sngfilter, wvlmin=wvlmin, wvlmax=wvlmax, gotdir=gotdir, file_abund=file_abund'
   return
endif

; Generates G(T,n) tables (200 x 200 pts) for the AIA response functions. 
; It assumes the chianti.ioneq CHIANTI file for ionization equilibrium
; values. 

; INPUT:
; sngfilter = (string) 'all' if all filters except UV (4500, 1700, 1600) are to be generated
;              'uv' if UV filters are to be generated (4500, 1700, 1600)
;              name of filter (ex: '171' for the AIA 171 filter):
;              EUV:
;              '304' -> AIA 304
;              '171' -> AIA 171
;              '193' -> AIA 193
;              '211' -> AIA 211
;              '335' -> AIA 335
;              '094' -> AIA 094
;              '131' -> AIA 131
;              UV:
;              '1600' -> AIA 1600
;              '1700' -> AIA 1700
;              '4500' -> AIA 4500
; wvlmin = (float) Minimum of desired wavelength range for line transition in Angstroms
; wvlmax = (float) Maximum of desired wavelength range for line transition in Angstroms
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

n_e_min = 1.e6
n_e_max = 1.e12
;n_e_min = 1.e8
;n_e_max = 1.e11

;steplg = 0.0015
steplg = 0.03
numn = alog10(n_e_max/n_e_min)/steplg

n_e_lg = dindgen(numn+1)/numn*alog10(n_e_max/n_e_min)+alog10(n_e_min)
n_e = 10.^(n_e_lg)

sterad_aia_pix=8.4d-12
watom = 0.

aia_resp = aia_get_response(/dn)
if sngfilter eq 'uv' then aia_resp_uv = aia_get_response(/dn,/uv)

if sngfilter eq 'all' or sngfilter eq '304' then begin
   openw,unit1,gotdir+'goft_table_aia304_'+nab+extname+'.dat',/get_lun & w0_1 = 304. & ion_1 = '304'
   printf,unit1,ion_1
   printf,unit1,w0_1
   printf,unit1,watom
   printf,unit1,numn,numt
   printf,unit1,alogt
endif
if sngfilter eq 'uv' or sngfilter eq '1600' then begin
   openw,unit2,gotdir+'goft_table_aia1600_'+nab+extname+'.dat',/get_lun & w0_2 = 1600. & ion_2 = '1600'
   printf,unit2,ion_2
   printf,unit2,w0_2
   printf,unit2,watom
   printf,unit2,numn,numt
   printf,unit2,alogt
endif
if sngfilter eq 'uv' or sngfilter eq '1700' then begin
   openw,unit3,gotdir+'goft_table_aia1700_'+nab+extname+'.dat',/get_lun & w0_3 = 1700. & ion_3 = '1700'
   printf,unit3,ion_3
   printf,unit3,w0_3
   printf,unit3,watom
   printf,unit3,numn,numt
   printf,unit3,alogt
endif
if sngfilter eq 'uv' or sngfilter eq '4500' then begin
   openw,unit4,gotdir+'goft_table_aia4500_'+nab+extname+'.dat',/get_lun & w0_4 = 4500. & ion_4 = '4500'
   printf,unit4,ion_4
   printf,unit4,w0_4
   printf,unit4,watom
   printf,unit4,numn,numt
   printf,unit4,alogt
endif
if sngfilter eq 'all' or sngfilter eq '171' then begin 
   openw,unit5,gotdir+'goft_table_aia171_'+nab+extname+'.dat',/get_lun & w0_5 = 171. & ion_5 = '171'
   printf,unit5,ion_5
   printf,unit5,w0_5
   printf,unit5,watom
   printf,unit5,numn,numt
   printf,unit5,alogt
endif
if sngfilter eq 'all' or sngfilter eq '193' then begin
   openw,unit6,gotdir+'goft_table_aia193_'+nab+extname+'.dat',/get_lun & w0_6 = 193. & ion_6 = '193'
   printf,unit6,ion_6
   printf,unit6,w0_6
   printf,unit6,watom
   printf,unit6,numn,numt
   printf,unit6,alogt
endif
if sngfilter eq 'all' or sngfilter eq '211' then begin
   openw,unit7,gotdir+'goft_table_aia211_'+nab+extname+'.dat',/get_lun & w0_7 = 211. & ion_7 = '211'
   printf,unit7,ion_7
   printf,unit7,w0_7
   printf,unit7,watom
   printf,unit7,numn,numt
   printf,unit7,alogt
endif
if sngfilter eq 'all' or sngfilter eq '335' then begin
   openw,unit8,gotdir+'goft_table_aia335_'+nab+extname+'.dat',/get_lun & w0_8 = 335. & ion_8 = '335'
   printf,unit8,ion_8
   printf,unit8,w0_8
   printf,unit8,watom
   printf,unit8,numn,numt
   printf,unit8,alogt
endif
if sngfilter eq 'all' or sngfilter eq '094' then begin 
   openw,unit9,gotdir+'goft_table_aia094_'+nab+extname+'.dat',/get_lun & w0_9 = 094. & ion_9 = '094'
   printf,unit9,ion_9
   printf,unit9,w0_9
   printf,unit9,watom
   printf,unit9,numn,numt
   printf,unit9,alogt
endif
if sngfilter eq 'all' or sngfilter eq '131' then begin 
   openw,unit10,gotdir+'goft_table_aia131_'+nab+extname+'.dat',/get_lun & w0_10 = 131. & ion_10 = '131'
   printf,unit10,ion_10
   printf,unit10,w0_10
   printf,unit10,watom
   printf,unit10,numn,numt
   printf,unit10,alogt
endif

for i=0,numn-1 do begin
   print,'doing density ',i,' of ',numn
   isothermal, wvlmin, wvlmax, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file

   if sngfilter eq 'all' or sngfilter eq '304' then eff_304=interpol(aia_resp.a304.ea,aia_resp.a304.wave,lambda,/spline)
   if sngfilter eq 'uv' or sngfilter eq '1600' then eff_1600=interpol(aia_resp_uv.a1600.ea,aia_resp.a1600.wave,lambda,/spline)
   if sngfilter eq 'uv' or sngfilter eq '1700' then eff_1700=interpol(aia_resp_uv.a1700.ea,aia_resp.a1700.wave,lambda,/spline)
   if sngfilter eq 'uv' or sngfilter eq '4500' then eff_4500=interpol(aia_resp_uv.a4500.ea,aia_resp.a4500.wave,lambda,/spline)
   if sngfilter eq 'all' or sngfilter eq '171' then eff_171=interpol(aia_resp.a171.ea,aia_resp.a171.wave,lambda,/spline)
   if sngfilter eq 'all' or sngfilter eq '193' then eff_193=interpol(aia_resp.a193.ea,aia_resp.a193.wave,lambda,/spline)
   if sngfilter eq 'all' or sngfilter eq '211' then eff_211=interpol(aia_resp.a211.ea,aia_resp.a211.wave,lambda,/spline)
   if sngfilter eq 'all' or sngfilter eq '335' then eff_335=interpol(aia_resp.a335.ea,aia_resp.a335.wave,lambda,/spline)
   if sngfilter eq 'all' or sngfilter eq '094' then eff_094=interpol(aia_resp.a94.ea,aia_resp.a94.wave,lambda,/spline)
   if sngfilter eq 'all' or sngfilter eq '131' then eff_131=interpol(aia_resp.a131.ea,aia_resp.a131.wave,lambda,/spline)
   
   if sngfilter eq 'all' or sngfilter eq '304' then begin sp_conv_304 = spectrum & sp_conv_304[*,*]=0. & endif
   if sngfilter eq 'uv' or sngfilter eq '1600' then begin sp_conv_1600 = spectrum & sp_conv_1600[*,*]=0. & endif
   if sngfilter eq 'uv' or sngfilter eq '1700' then begin sp_conv_1700 = spectrum & sp_conv_1700[*,*]=0. & endif
   if sngfilter eq 'uv' or sngfilter eq '4500' then begin sp_conv_4500 = spectrum & sp_conv_4500[*,*]=0. & endif
   if sngfilter eq 'all' or sngfilter eq '171' then begin sp_conv_171 = spectrum & sp_conv_171[*,*]=0. & endif
   if sngfilter eq 'all' or sngfilter eq '193' then begin sp_conv_193 = spectrum & sp_conv_193[*,*]=0. & endif
   if sngfilter eq 'all' or sngfilter eq '211' then begin sp_conv_211 = spectrum & sp_conv_211[*,*]=0. & endif
   if sngfilter eq 'all' or sngfilter eq '335' then begin sp_conv_335 = spectrum & sp_conv_335[*,*]=0. & endif
   if sngfilter eq 'all' or sngfilter eq '094' then begin sp_conv_094 = spectrum & sp_conv_094[*,*]=0. & endif
   if sngfilter eq 'all' or sngfilter eq '131' then begin sp_conv_131 = spectrum & sp_conv_131[*,*]=0. & endif
   
   for j=0,numt-1 do begin
      if sngfilter eq 'all' or sngfilter eq '304' then sp_conv_304[*,j]=sterad_aia_pix*spectrum[*,j]*eff_304
      if sngfilter eq 'uv' or sngfilter eq '1600' then sp_conv_1600[*,j]=sterad_aia_pix*spectrum[*,j]*eff_1600
      if sngfilter eq 'uv' or sngfilter eq '1700' then sp_conv_1700[*,j]=sterad_aia_pix*spectrum[*,j]*eff_1700
      if sngfilter eq 'uv' or sngfilter eq '4500' then sp_conv_4500[*,j]=sterad_aia_pix*spectrum[*,j]*eff_4500
      if sngfilter eq 'all' or sngfilter eq '171' then sp_conv_171[*,j]=sterad_aia_pix*spectrum[*,j]*eff_171
      if sngfilter eq 'all' or sngfilter eq '193' then sp_conv_193[*,j]=sterad_aia_pix*spectrum[*,j]*eff_193
      if sngfilter eq 'all' or sngfilter eq '211' then sp_conv_211[*,j]=sterad_aia_pix*spectrum[*,j]*eff_211
      if sngfilter eq 'all' or sngfilter eq '335' then sp_conv_335[*,j]=sterad_aia_pix*spectrum[*,j]*eff_335
      if sngfilter eq 'all' or sngfilter eq '094' then sp_conv_094[*,j]=sterad_aia_pix*spectrum[*,j]*eff_094
      if sngfilter eq 'all' or sngfilter eq '131' then sp_conv_131[*,j]=sterad_aia_pix*spectrum[*,j]*eff_131
   endfor
   if sngfilter eq 'all' or sngfilter eq '304' then   resp_304=total(sp_conv_304,1)
   if sngfilter eq 'uv' or sngfilter eq '1600' then   resp_1600=total(sp_conv_1600,1)
   if sngfilter eq 'uv' or sngfilter eq '1700' then   resp_1700=total(sp_conv_1700,1)
   if sngfilter eq 'uv' or sngfilter eq '4500' then   resp_4500=total(sp_conv_4500,1)
   if sngfilter eq 'all' or sngfilter eq '171' then   resp_171=total(sp_conv_171,1)
   if sngfilter eq 'all' or sngfilter eq '193' then   resp_193=total(sp_conv_193,1)
   if sngfilter eq 'all' or sngfilter eq '211' then   resp_211=total(sp_conv_211,1)
   if sngfilter eq 'all' or sngfilter eq '335' then   resp_335=total(sp_conv_335,1)
   if sngfilter eq 'all' or sngfilter eq '094' then   resp_094=total(sp_conv_094,1)
   if sngfilter eq 'all' or sngfilter eq '131' then   resp_131=total(sp_conv_131,1)

   if sngfilter eq 'all' or sngfilter eq '304' then begin
      printf,unit1,n_e_lg[i]
      printf,unit1,resp_304
   endif
   if sngfilter eq 'uv' or sngfilter eq '1600' then begin
      printf,unit2,n_e_lg[i]
      printf,unit2,resp_1600
   endif
   if sngfilter eq 'uv' or sngfilter eq '1700' then begin
      printf,unit3,n_e_lg[i]
      printf,unit3,resp_1700
   endif
   if sngfilter eq 'uv' or sngfilter eq '4500' then begin
      printf,unit4,n_e_lg[i]
      printf,unit4,resp_4500
   endif
   if sngfilter eq 'all' or sngfilter eq '171' then begin
      printf,unit5,n_e_lg[i]
      printf,unit5,resp_171
   endif
   if sngfilter eq 'all' or sngfilter eq '193' then begin
      printf,unit6,n_e_lg[i]
      printf,unit6,resp_193
   endif
   if sngfilter eq 'all' or sngfilter eq '211' then begin
      printf,unit7,n_e_lg[i]
      printf,unit7,resp_211
   endif
   if sngfilter eq 'all' or sngfilter eq '335' then begin
      printf,unit8,n_e_lg[i]
      printf,unit8,resp_335
   endif
   if sngfilter eq 'all' or sngfilter eq '094' then begin
      printf,unit9,n_e_lg[i]
      printf,unit9,resp_094
   endif
   if sngfilter eq 'all' or sngfilter eq '131' then begin
      printf,unit10,n_e_lg[i]
      printf,unit10,resp_131
   endif
   print,string(13b)+' % finished: ',float(i)*100./(numn-1),format='(a,f4.0,$)'
endfor

if sngfilter eq 'all' or sngfilter eq '304' then free_lun,unit1
if sngfilter eq 'uv' or sngfilter eq '1600' then free_lun,unit2
if sngfilter eq 'uv' or sngfilter eq '1700' then free_lun,unit3
if sngfilter eq 'uv' or sngfilter eq '4500' then free_lun,unit4
if sngfilter eq 'all' or sngfilter eq '171' then free_lun,unit5
if sngfilter eq 'all' or sngfilter eq '193' then free_lun,unit6
if sngfilter eq 'all' or sngfilter eq '211' then free_lun,unit7
if sngfilter eq 'all' or sngfilter eq '335' then free_lun,unit8
if sngfilter eq 'all' or sngfilter eq '094' then free_lun,unit9
if sngfilter eq 'all' or sngfilter eq '131' then free_lun,unit10

end
