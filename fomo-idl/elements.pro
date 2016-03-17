
pro elements,w0=w0,ion=ion,logTm=logTm,enum=enum,inum=inum,ind=ind,vers=vers

if keyword_set(w0) eq 0 then begin
   print,'elements,w0=w0,ion=ion,logTm=logTm,enum=enum,inum=inum,ind=ind,vers=vers'
   return
endif

; Reads info for definition of line transitions

; INPUT:
; w0: wavelength center of line in Angstroms (float)
; ion: Ion of interest (string, eg:'fe_9')
; vers: Chianti version. Only version 6 or 7. If higher then the index
; number needs to be calculated as explained below. 

; OUTPUT:
; logTm: Temperature of maximum emission 
; (can be calculated with:
; integral_calc,enum,inum,ind=ind,outstr=outstr) ; assumes log(dens)=10
; logTm = outstr.tmem
; enum: nuclear charge of element
; inum: Ionization number
; ind: index for wavelength of transition (given by emiss_calc &
; emiss_select function)
; ex: emiss = emiss_calc(enum,inum)
;     n_elt = n_elements(emiss.lambda)
;     index = lindgen(n_elt)
;     wavels = emiss(index).lambda
;     wrange = [w0-dw,w0+dw] ;dw is a sufficiently small wavelength
;     so that the wavelength interval only contains the wanted transition
;     ind = index[0]
;     another way is like this:
;     cal_emiss = emiss_select(emiss,wrange=[170,172],sel_ind=ind)
 
case ion of
; for Fe IX 171.073  
   'fe_9': begin
      case w0 of
         171.073: begin
            logTm = 5.93         ; max formation Temperature of line
            enum = 26           ; number for element
            inum = 9            ; number of ion
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 151 
               if vers eq 7 then ind = 5870 ; index for wavelength of transition
            endif
         end
      endcase
   end
; for Fe XII 193.509
   'fe_12': begin
      case w0 of 
         193.509: begin
            logTm = 6.19
            enum = 26
            inum = 12
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 844
               if vers eq 7 then ind = 23050
            endif
         end
         195.12: begin ;3s2 3p3 4S3/2 - 3s2 3p2 (3P)  3d 4P5/2
            logTm = 6.19
            enum = 26
            inum = 12
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 873
               if vers eq 7 then ind = 23391
            endif
         end
         186.88: begin
            logTm = 6.19
            enum = 26
            inum = 12
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 750
               if vers eq 7 then ind = 21747
            endif
         end
      endcase
   end
; for Fe XIII 202 and 203:
   'fe_13': begin
      case w0 of
         202.044: begin
            logTm = 6.25
            enum = 26
            inum = 13
            if n_elements(vers) ne 0 then begin
               if vers eq 7 then ind = 39127
            endif
         end
         203.828: begin
            logTm = 6.26
            enum = 26
            inum = 13
            if n_elements(vers) ne 0 then begin
               if vers eq 7 then ind = 39486
            endif
         end
      endcase
   end
   'fe_15': begin
      case w0 of
          284.16: begin
	    enum = 26
	    inum = 15
            if n_elements(vers) ne 0 then begin
               if vers eq 7 then begin
                  ind = 5609
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
	    logTm = 6.34
          end
      endcase
   end
; Si VII 197.7684
   'si_7': begin
      case w0 of
         197.768: begin
            enum = 14
            inum = 7
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 591 
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            logTm = 5.8
         end
      endcase
   end
; He II 303.781 & 303.786
   'he_2':begin
      case w0 of
         303.781: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 6 
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 2
            inum = 2
            logTm = 4.92
         end
         303.786: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 7
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 2
            inum = 2
            logTm = 4.92
         end
      endcase
   end
; C IV 1548.189 & 1550.775
   'c_4': begin
      case w0 of
         1548.19: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 16
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 6
            inum = 4
            logTm = 5.04   
         end
         1550.78: begin ;1550.775
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 17
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 6
            inum = 4
            logTm = 5.04        
         end
      endcase
   end
; Ne VIII 770.4103 
   'ne_8': begin
      case w0 of
         770.41: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 86
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse  
            endif
            enum = 10
            inum = 8
            logTm = 5.8
         end
      endcase
   end
; O IV  1399.78 & 1401.16
   'o_4': begin
      case w0 of
         1399.77: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 512
               if vers eq 7 then ind = 716
            endif
            enum = 8
            inum = 4
            logTm = 5.17
         end
         1401.16: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 513
               if vers eq 7 then ind = 717
            endif
            enum = 8
            inum = 4
            logTm = 5.17
         end
         1404.78: begin;1404.81?
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then ind = 514 
               if vers eq 7 then ind = 718
            endif
            enum = 8
            inum = 4
            logTm = 5.17
         end
      endcase
   end
; Si IV 1393.757 & 1402.772
   'si_4': begin
      case w0 of
         1393.76: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 26
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 14
            inum = 4
            logTm = 4.89
         end
         1402.77: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 27
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 14
            inum = 4
            logTm = 4.89
         end
      endcase
   end
; C III 977.02 & 1174.933
   'c_3': begin
      case w0 of
         977.02: begin;2s2 1S0 - 2s.2p 1P1
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 21
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 6
            inum = 3
            logTm = 4.93
         end
         1174.93: begin;2s.2p 3P1 - 2p2 3P2
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 24
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 6
            inum = 3
            logTm = 4.94
         end
      endcase
   end
; C II 1334.535 & 1335.71
   'c_2': begin
      case w0 of
         1334.53: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 then begin
                  ind = 18
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 6
            inum = 2
            logTm = 4.74
         end
         1335.71: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 56
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 6
            inum = 2
            logTm = 4.62
         end
      endcase
   end
; Mg II h & k: 2796.3521 & 2803.531
   'mg_2': begin
      case w0 of
         2796.35: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 11
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 12
            inum = 2
            logTm = 4.22
         end
         2803.531: begin
            if n_elements(vers) ne 0 then begin
               if vers eq 6 or vers eq 7 then begin
                  ind = 14
               endif else begin
                  print,'check index corresponding to this Chianti version'
                  stop
               endelse
            endif
            enum = 12
            inum = 2
            logTm = 4.22
         end
      endcase
   end
endcase

end
