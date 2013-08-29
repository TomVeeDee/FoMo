
pro lookup_goft, ion=ion, w0=w0, gotdir=gotdir,n_e_lg=n_e_lg, logt=logt, goft_mat=goft_mat, watom=watom, filenm=filenm,file_abund=file_abund

if keyword_set(ion) eq 0 then begin
   print,'Check input directories'
   print,'lookup_goft, ion=ion, w0=w0, n_e_lg=n_e_lg, logt=logt, goft_mat=goft_mat, watom=watom'
   return
endif

; returns got_mat matrix reading dat files created by table_goft.pro
; INPUT:
; ion = spectral line : e.g. 'fe_9' for Fe IX
; w0: wavelength center of line
; OUTPUT:
; n_e_lg: number density array (logarithm) of table
; logt: temperature array (logarithm) of table
; goft_mat: G(T,n) function evaluated at (n_e, t)

; INPUT DIRECTORY:
;gotdir = '/users/cpa/pantolin/Modeling/FoMo/chiantitables/'
w0nm = string(round(w0),format='(i4.4)')
if keyword_set(file_abund) then begin
   if file_abund eq 'coronal' then nab = '_abco'
   if file_abund eq 'photospheric' then nab = '_abph'
endif else begin
   nab = '_abco'
endelse

if keyword_set(filenm) then filegot = 'goft_table_'+filenm+ion+nab+'.dat' else filegot = 'goft_table_'+ion+'_'+w0nm+nab+'.dat'

;if ion eq 'fe_9' then filegot = 'goft_table_f2rt_171.dat'
;if ion eq 'fe_12' then filegot = 'goft_table_frt_193.dat'

if file_test(gotdir+filegot) eq 0 then filename = dialog_pickfile(filter='*.dat',/read) else filename = gotdir+filegot

print,'G(n,T) table: '+filegot
;  n_e_min = 1.e8
;  n_e_max = 1.e11 ; for goft_table_frt_193.dat or goft_table_f2rt_171.dat choose n_e_max = 1.e10
;  steplg = 0.001
;  num = alog10(n_e_max/n_e_min)/steplg

  watom = 0
  openr,unit,filename,/get_lun
  readf,unit,ion
  readf,unit,w00
  readf,unit,watom
  readf,unit,numn,numt
  logt = fltarr(numt)
  readf,unit,logt
  goft_mat = fltarr(numn,numt)
  n_e_lg = dblarr(numn)
  g_t = fltarr(numt)
  i = 0

  while(eof(unit) eq 0) do begin
     readf,unit,n_e_0
     readf,unit,g_t
     goft_mat[i,*] = g_t
     n_e_lg[i] = n_e_0
     i = i+1
  endwhile
  free_lun,unit
end
