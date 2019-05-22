pro make_eitresponse, sngfilter=sngfilter,wvlmin=wvlmin,wvlmax=wvlmax,gotdir=gotdir,file_abund=file_abund,extname=extname,extro=extro,silent=silent

if ~keyword_set(sngfilter) then begin
   print,'make_eitresponse, sngfilter=sngfilter,wvlmin=wvlmin,wvlmax=wvlmax,gotdir=gotdir,file_abund=file_abund,extname=extname,extro=extro,silent=silent'
   return
endif

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
; extro: (string) for G(T,n) tables with extended density range [6,12]
;        in log. Default is [8,12] for EIT 304 and [8,11] for the rest

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

if keyword_set(extro) then begin
   n_e_min = 1.e6
   n_e_max_sml = 1.e12
   n_e_max_big = 1.e12
   extname = extname+'_extro'
endif else begin
   n_e_min = 1.e8
   n_e_max_sml = 1.e11
   n_e_max_big = 1.e12
endelse

steplg = 0.015
numn_sml = round(alog10(n_e_max_sml/n_e_min)/steplg)
numn_big = round(alog10(n_e_max_big/n_e_min)/steplg)

n_e_lg_sml = dindgen(numn_sml+1)/numn_sml*alog10(n_e_max_sml/n_e_min)+alog10(n_e_min)
n_e_lg_big = dindgen(numn_big+1)/numn_big*alog10(n_e_max_big/n_e_min)+alog10(n_e_min)
n_e_sml = 10.^(n_e_lg_sml)   
n_e_big = 10.^(n_e_lg_big)  

;sterad_eit_pix=(1.275e-5)^2 ;;;  EIT pixel - 2.629arcsec; sterad_eit_pix=(1.275e-5)^2=(2.629arcsec)^2
sterad_eit_pix=1.0
watom = 0.
unitstring = 'cm^5 DN s^{-1} sr^{-1}'
openr,unitversion, !xuvtop+'/VERSION',/get_lun
str=''
readf,unitversion, str
free_lun,unitversion
vchianti = 'CHIANTI'+str

if sngfilter[0] eq 'all' then eitarr = ['171','195','284','304'] else eitarr = sngfilter
neitar = n_elements(eitarr)
NN = 1000
wmin171=160.0 & wmax171=250.0 & lambda171=findgen(NN)/(NN-1)*(wmax171-wmin171)+wmin171
wmin195=160.0 & wmax195=250.0 & lambda195=findgen(NN)/(NN-1)*(wmax195-wmin195)+wmin195
wmin284=200.0 & wmax284=330.0 & lambda284=findgen(NN)/(NN-1)*(wmax284-wmin284)+wmin284
wmin304=200.0 & wmax304=350.0 & lambda304=findgen(NN)/(NN-1)*(wmax304-wmin304)+wmin304

for k=0,neitar-1 do begin
   filt = eitarr[k]
   print,'Constructing EIT - '+filt+' G(T) table'
   ex = 'eit_resp_'+filt+' = eit_parms(lambda'+filt+',fix(filt),units=units)'
   void = execute(ex)
   openw,unit,gotdir+'goft_table_eit'+filt+'_'+nab+extname+'.dat',/get_lun & w0 = float(filt) & ion = filt
   ex = 'wvlmin = wmin'+filt
   void = execute(ex)
   ex = 'wvlmax = wmax'+filt
   void = execute(ex)
   if filt eq '304' then begin
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
   printf,unit,wvlmin,wvlmax
   printf,unit,unitstring
   printf,unit,vchianti
   printf,unit,numn,numt
   printf,unit,alogt

   for i=0,numn-1 do begin
      print,'EIT '+filt+' - doing density number ',i,' of ',numn
      isothermal, wvlmin, wvlmax, 0.1, temp, lambda, spectrum, list_wvl,list_ident, edensity=n_e[i] ,/photons,/cont,/all,ioneq_name=ioneq_name,abund_name=abund_file,noverbose=silent
      ex = 'eff = interpol(eit_resp_'+filt+',lambda'+filt+',lambda)'
      void = execute(ex)
      sp_conv = spectrum & sp_conv[*,*]=0. 
      for j=0,numt-1 do sp_conv[*,j]=fnhne*sterad_eit_pix*spectrum[*,j]*eff
      resp = total(sp_conv,1)
      printf,unit,n_e_lg[i]
      printf,unit,resp
   endfor
   free_lun,unit
endfor

end
