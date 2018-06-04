
pro lookup_goft,ion=ion,w0=w0, gotdir=gotdir,n_e_lg=n_e_lg, logt=logt, goft_mat=goft_mat, watom=watom,filenm=filenm,channel=channel,nab=nab,silent=silent

if keyword_set(ion) eq 0 then begin
   print,'Check input directories'
   print,'lookup_goft, ion=ion, w0=w0, gotdir=gotdir,n_e_lg=n_e_lg, logt=logt, goft_mat=goft_mat, watom=watom, filenm=filenm,channel=channel,nab=nab,silent=silent'
   return
endif

; Returns G(T,n) values reading dat files created by table_goft.pro

; INPUT:
; ion: (string) acronym of the ion (not needed for channels like SDO/AIA filter)
; w0: (float) wavelength of line center (filter number for SDO/AIA filter)
; gotdir: (string) directory path to where the .dat G(T,n) file is.
; filenm: (string) name helping defining the .dat file

; OUTPUT:
; n_e_lg: number density array (logarithm) where G(T,n) is defined
; logt: temperature array (logarithm) where G(T,n) is defined
; goft_mat: G(T,n) function evaluated at (n_e, t)
; watom = get_atomic_weight(enum), where enum is the nuclear charge of element

if keyword_set(w0) then begin
   if w0 lt 1.e3 and keyword_set(channel) then w0nm = string(round(w0),format='(i3.3)')
   if w0 lt 1.e3 and ~keyword_set(channel) then w0nm = string(round(w0),format='(i4.4)')
   if w0 gt 1.e3 and w0 lt 1.e4 then w0nm = string(round(w0),format='(i4.4)')
   if w0 gt 1.e4 then w0nm = string(round(w0),format='(i5)')
endif

; by default looks for tables with coronal abundances:
if ~keyword_set(nab) then begin
   print,'Abundance package not specified. Adopting coronal abundances.'
   nab = '_abco'
endif

if keyword_set(filenm) then filegot = 'goft_table_'+filenm+w0nm+nab+'.dat' 
if keyword_set(channel) then filegot = 'goft_table_'+channel+w0nm+nab+'.dat' 
if (~keyword_set(filenm) and ~keyword_set(channel)) then filegot = 'goft_table_'+ion+'_'+w0nm+nab+'.dat'

if file_test(gotdir+filegot) eq 0 then begin
   print,'G(T,n) table could not be found'
   filename = dialog_pickfile(filter='*.dat',/read) 
endif else begin 
   filename = gotdir+filegot
endelse
if ~keyword_set(silent) then print,'G(n,T) table: '+filegot

  watom = 0 & lv1 = 0 & lv2 = 0

  openr,unit,filename,/get_lun
  readf,unit,ion
  readf,unit,cw0
  readf,unit,watom
  readf,unit,lv1
  readf,unit,lv2
  readf,unit,numn,numt
  logt = fltarr(numt)
  readf,unit,logt
  goft_mat = fltarr(numn,numt)
  n_e_lg = dblarr(numn)
  g_t = fltarr(numt)
  i = 0

  if w0 ne cw0 then begin
     print,'Warning: input wavelength (Angs.):',w0
     print,'G(T) wavelength (Angs.): ',cw0
  endif
  if ~keyword_set(quiet) then begin
     if keyword_set(channel) then begin
        print,'Channel: ',ion
        print,'Representative wavelength (Angs.): ',cw0
        wvrng = strcompress('['+string(lv1)+','+string(lv2)+']',/remove_all)
        print,'Wavelength range of channel (FHWHM in Angs.): '+wvrng
     endif
     if ~keyword_set(channel) then begin
        print,'Ion of line transition: ',ion
        print,'Wavelength of line transition (Angs.): ',cw0
        print,'Atomic weight: ',watom
        print,'Index of lower level of line transition: ',lv1
        print,'Index of upper level of line transition: ',lv2
     endif
  endif

  while(eof(unit) eq 0) do begin
     readf,unit,n_e_0
     readf,unit,g_t
     goft_mat[i,*] = g_t
     n_e_lg[i] = n_e_0
     i = i+1
  endwhile
  free_lun,unit
end
