
pro elements,ion=ion,w0=w0,cw0=cw0,enum=enum,inum=inum,ind=ind,lv1=lv1,lv2=lv2,silent=silent

  if ~keyword_set(ion) then begin
     print,'ion=ion,w0=w0,cw0=cw0,enum=enum,inum=inum,ind=ind,lv1=lv1,lv2=lv2,silent=silent'
     return
  endif
; Reads info for definition of line transitions

; INPUT:
; w0: wavelength center of line in Angstroms (float)
; ion: Ion of interest (string, eg:'fe_9')

; OUTPUT:
; w0: wavelength (in Amstrogms) of line transition
; cw0: corresponding closest wavelength to w0 in the CHIANTI database
; enum: nuclear charge of element
; inum: Ionization number
; ind: index for wavelength of transition (given by emiss_calc &
; emiss_select function)
; lv1 = index of lower level of line transition (same
; nomenclature as Chianti's read_wgfa2)
; lv2 = index of upper level of line transition (same
; nomenclature as Chianti's read_wgfa2)
  
; examples:
; Fe IX 171.073
; ion = 'fe_9'
; w0 = 171.07
; cw0 = 171.073
; enum = 26 
; inum = 9  
; lv1 = 1
; lv2 = 13

        convertname,ion,enum,inum ; calculates ion number enum and ionization number inum
        emiss = emiss_calc(enum,inum,/no_calc,quiet=silent) ; calculates all possible lines
        wavels = emiss.lambda
        lev1 = emiss.level1
        lev2 = emiss.level2
        ind = ([min(abs(wavels - w0)),!c])[1]
        cw0 = wavels[ind]
        lv1 = lev1[ind]
        lv2 = lev2[ind]        
end
