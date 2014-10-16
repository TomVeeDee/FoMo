

pro goft_table, w0=w0, ion=ion,gotdir=gotdir,vers=vers

if keyword_set(w0) eq 0 or keyword_set(ion) eq 0 then begin
   print,'goft_table, w0=w0, ion=ion,gotdir=gotdir,vers=vers'
   return
endif

; Generates a G(T,n) table (200 pts in temperature, 3000 points in density) for a given line transition.
; A coronal abundance is included by default ('sun_coronal.abund'). Such abundance can be replaced
; or eliminated in the routine lineongrid_goft_tab.pro

; INPUT:
; w0 = (float) line center wavelength of line transition in Angstroms
; ion = (string) acronym of ion (see elements.pro)
; gotdir = (string) directory path where to save the generated table
;          (don't forget '/' at end of path) 
; vers = (int) Chianti version (6 or 7)

; CALLS:
; elements, get_atomic_weight, g_of_t

if ~keyword_set(vers) then begin
   vers = 7
   print,'Assuming Chianti version 7'
endif

abund_file = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal.abund')
nab = 'abco'

; Given wavelength center of line (w0) and ion of interest (ion)
; produces the G(T,n) table for that wavelength transition

  ; determine wavelength transition:
  elements, w0=w0, ion=ion, logTm=logTm, enum=enum, inum=inum, ind=ind, vers=vers

  ; get atomic weight:
  watom = get_atomic_weight(enum)

  n_e_min = 1.e8
  n_e_max = 1.e11

steplg = 0.001

numn = alog10(n_e_max/n_e_min)/steplg

n_e_lg = dindgen(numn+1)/numn*alog10(n_e_max/n_e_min)+alog10(n_e_min)

w0nm = string(round(w0),format='(i4.4)')
openw,unit,gotdir+'goft_table_'+ion+'_'+w0nm+'_'+nab+'.dat',/get_lun
printf,unit,ion
printf,unit,w0
printf,unit,watom

alogt0 = findgen(101)/20+4. ; same range as that defined by Chianti
wte = 0.75
tmin = (logTm-wte)>4.           ;define wanted range
tmax = (logTm+wte)<7.
pts = where(alogt0 ge tmin and alogt0 le tmax,npts)
alogt1 = alogt0[pts]
numt = 200
alogt2 = findgen(numt)/(numt-1)*2*wte+tmin

printf,unit,numn,numt
printf,unit,alogt2
for i=0,numn do begin
   goft0=g_of_t(enum,inum,dens=n_e_lg[i],ioneq_file=concat_dir(concat_dir(!xuvtop,'ioneq'),'chianti.ioneq'),abund_file=abund_file,index=ind,/quiet)

   goft1 = goft0[pts]
   ion_interp,alogt1,goft1,alogt2,goft2

   printf,unit,n_e_lg[i]
   printf,unit,goft2
endfor
free_lun,unit

end
