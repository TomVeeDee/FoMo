
pro elements,ion=ion,w0=w0,enum=enum,inum=inum,ind=ind

  if ~keyword_set(ion) then begin
     print,'elements,ion=ion,w0=w0,enum=enum,inum=inum,ind=ind'
     return
  endif
; Reads info for definition of line transitions

; INPUT:
; w0: wavelength center of line in Angstroms (float)
; ion: Ion of interest (string, eg:'fe_9')

; OUTPUT:
; w0: w0 is updated to the wavelength as found by CHIANTI
; enum: nuclear charge of element
; inum: Ionization number
; ind: index for wavelength of transition (given by emiss_calc &
; emiss_select function)

; examples:
; Fe IX 171.073
; ion = 'fe_9'
; w0 = 171.073
; enum = 26 (number for element)
; inum = 9  (number of ion)

; for Fe XII 193.509
;  ion = 'fe_12'
;  w0 = 193.509
;  logTm = 6.19
;  enum = 26
;  inum = 12

;  ion = 'fe_12'
;  w0 = 195.12 (3s2 3p3 4S3/2 - 3s2 3p2 (3P)  3d 4P5/2)
;  logTm = 6.19
;  enum = 26
;  inum = 12

;  ion = 'fe_12'
;  w0 = 186.88
;  logTm = 6.19
;  enum = 26
;  inum = 12

; for Fe XIII 202 and 203:
;  ion = 'fe_13'
;  w0 = 202.044
;  logTm = 6.25
;  enum = 26
;  inum = 13

;  ion = 'fe_13'
;  w0 = 203.828
;  logTm = 6.26
;  enum = 26
;  inum = 13

;  ion = 'fe_15'
;  w0 = 284.16
;  enum = 26
;  inum = 15
;  logTm = 6.34

; Si VII 197.7684
;  ion = 'si_7'
;  w0 = 197.768
;  enum = 14
;  inum = 7
;  logTm = 5.8

; He II 303.781 & 303.786
;  ion = 'he_2'
;  w0 = 303.781
;  enum = 2
;  inum = 2
;  logTm = 4.92

;  ion = 'he_2'
;  w0 = 303.786
;  enum = 2
;  inum = 2
;  logTm = 4.92

; C IV 1548.189 & 1550.775
;  ion = 'c_4'
;  w0 = 1548.19
;  enum = 6
;  inum = 4
;  logTm = 5.04   

;  ion = 'c_4'
;  w0 = 1550.78
;  enum = 6
;  inum = 4
;  logTm = 5.04        

; Ne VIII 770.4103 
;  ion = 'ne_8'
;  w0 = 770.41
;  enum = 10
;  inum = 8
;  logTm = 5.8

; O IV  1399.78 & 1401.16
; ion = 'o_4'
; w0 = 1399.77
; enum = 8
; inum = 4
; logTm = 5.17

;  ion = 'o_4'
;  w0 = 1401.16
;  enum = 8
;  inum = 4
;  logTm = 5.17

;  ion = 'o_4'
;  w0 = 1404.78
;  enum = 8
;  inum = 4
;  logTm = 5.17

; Si IV 1393.757 & 1402.772
;  ion = 'si_4'
;  w0 = 1393.76
;  enum = 14
;  inum = 4
;  logTm = 4.89

;  ion = 'si_4'
;  w0 = 1402.77
;  enum = 14
;  inum = 4
;  logTm = 4.89

; C III 977.02 & 1174.933
;  ion = 'c_3'
;  w0 = 977.02  ;2s2 1S0 - 2s.2p 1P1
;  enum = 6
;  inum = 3
;  logTm = 4.93

;  ion = 'c_3'
;  w0 = 1174.93   ;2s.2p 3P1 - 2p2 3P2
;  enum = 6
;  inum = 3
;  logTm = 4.94

; C II 1334.535 & 1335.71
;  ion = 'c_2'
;  w0 = 1334.53
;  enum = 6
;  inum = 2
;  logTm = 4.74

;  ion = 'c_2'
;  w0 = 1335.71
;  enum = 6
;  inum = 2
;  logTm = 4.62

; Mg II h & k: 2796.3521 & 2803.531
;  ion = 'mg_2'
;  w0 = 2796.35
;  enum = 12
;  inum = 2
;  logTm = 4.22

;  ion = 'mg_2'
;  w0 = 2803.531
;  enum = 12
;  inum = 2
;  logTm = 4.22

  
        convertname,ion,enum,inum ; calculates ion number enum and ionization number inum
        emiss=emiss_calc(enum,inum,/no_calc,/quiet) ; calculates all possible lines
        wavels=emiss.lambda
        ind = ([min(abs(wavels - w0)),!c])[1]
	w0=wavels[ind]
        
end
