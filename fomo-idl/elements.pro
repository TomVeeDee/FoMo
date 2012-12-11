
pro elements, w0=w0, ion=ion, logTm=logTm, enum=enum, inum=inum, ind=ind

if keyword_set(w0) eq 0 then begin
   print,'elements, w0=w0, ion=ion, logTm=logTm, enum=enum, inum=inum, ind=ind'
   return
endif

; INPUT:
; w0: wavelength center of line (float)
; ion: Ion of interest (string, eg:'fe_9')
; OUTPUT:
; logTm: maximum formation Temperature of line
; enum: nuclear charge of element
; inum: Ionization number
; ind: index for wavelength of transition (given by emiss_calc function)

case ion of
; for Fe IX 171.073  
   'fe_9': begin
      case w0 of
         171.073: begin
            logTm = 5.9         ; max formation Temperature of line
            enum = 26           ; number for element
            inum = 9            ; number of ion
            ind = 151           ; index for wavelength of transition
         end
      endcase
   end
; for Fe XII 193.509
   'fe_12': begin
      case w0 of 
         193.509: begin
            logTm = 6.2
            enum = 26
            inum = 12
            ind = 844
         end
      endcase
   end
; Si VII 197.7684
   'si_7': begin
      case w0 of
         197.7684: begin
            enum = 14
            inum = 7
            ind = 591
            logTm = 5.9
         end
      endcase
   end
; He II 303.781 & 303.786
   'he_2':begin
      case w0 of
         303.781: begin
            ind = 6 
            enum = 2
            inum = 2
            logTm = 5.0
         end
         303.786: begin
            ind = 7
            enum = 2
            inum = 2
            logTm = 5.0
         end
      endcase
   end
; C IV 1548.189 & 1550.775
   'c_4': begin
      case w0 of
         1548.189: begin
            ind = 16
            enum = 6
            inum = 4
            logTm = 5.0         
         end
         1550.775: begin
            ind = 17
            enum = 6
            inum = 4
            logTm = 5.0         
         end
      endcase
   end
; Ne VIII 770.4103 
   'ne_8': begin
      case w0 of
         770.4103: begin
            ind = 86
            enum = 10
            inum = 8
            logTm = 5.9
         end
      endcase
   end
; O IV  1399.78 & 1401.16
   'o_4': begin
      case w0 of
         1399.78: begin
            ind = 512
            enum = 8
            inum = 4
            logTm = 5.1
         end
         1401.16: begin
            ind = 513
            enum = 8
            inum = 4
            logTm = 5.1
         end
      endcase
   end
; Si IV 1393.757 & 1402.772
   'si_4': begin
      case w0 of
         1393.757: begin
            ind = 26
            enum = 14
            inum = 4
            logTm = 4.8
         end
         1402.772: begin
            ind = 27
            enum = 14
            inum = 4
            logTm = 4.8
         end
      endcase
   end
; C II 1334.535 & 1335.71
   'c_2': begin
      case w0 of
         1334.535: begin
            ind = 18
            enum = 6
            inum = 2
            logTm = 4.0
         end
         1335.71: begin
            ind = 20
            enum = 6
            inum = 2
            logTm = 4.0
         end
      endcase
   end
; Mg II h & k: 2796.3521 & 2803.531
   'mg_2': begin
      case w0 of
         2796.3521: begin
            ind = 11
            enum = 12
            inum = 2
            logTm = 4.0
         end
         2803.531: begin
            ind = 14
            enum = 12
            inum = 2
            logTm = 4.0
         end
      endcase
   end
endcase

end
