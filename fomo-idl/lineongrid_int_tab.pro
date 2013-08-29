pro lineongrid_int_tab, rh_s, te_s, dirgot=dirgot,wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, logt=logt,watom=watom,intens=intens,conv=conv

if n_params(0) lt 1 then begin
   print,'Check input directories'
   print,'lineongrid_goft_tab, rh_s, te_s, wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion, w0=w0, intens=intens'
   return
endif

; Calculates the emission at each point of a given numerical box by
; reading tabulated G(T) values produced by function goft_table.pro
; INPUT:
; rh_j = 2D or 3D array of mass density in kg/m^3, normalized by 1.e10
; te_j = 2D or 3D array of temperature in K, normalized by
; (protonmass/(2*kboltz)*1.e10)
; OPTIONAL:
; ion = ion for which to calculate emissivities (default = fe_9)
; w0 = wavelength center of line transition (default = 171.073)
; minwave = infimum of wavelength range (default = 171.0)
; maxwave = suppremum of wavelength range (default = 171.14)
; OUTPUT:
; wave = wavelength array set to nwave pts, containing line transition	
; w0 = wavelength center
; emission_goft = array of emission (x, y, (z), lambda)

; set ionization,abundance and DEM packages for CHIANTI:
;ioneq_name = '/users/cpa/tomvd/ssw/packages/chianti/dbase/ioneq/chianti.ioneq'
;abund_name = '/users/cpa/tomvd/ssw/packages/chianti/dbase/abundance/sun_coronal.abund'
;dem_name = '/users/cpa/tomvd/ssw/packages/chianti/dbase/dem/active_region.dem'
ioneq_name= !xuvtop+'/ioneq/chianti.ioneq'
abund_name=!xuvtop+'/abundance/sun_photospheric.abund'
abund_name=!xuvtop+'/abundance/sun_coronal.abund'
;dem_name=!xuvtop+'/dem/active_region.dem'

proton=1.67262158*10^(-27.)
kboltz = 1.380658*10^(-23.)
gamma = 5./3.
mu = 1.27

normro = 1.e10
;   const = 3.5991482e+13 ; = T[mid]*rho[mid]^(-gamma+1)
;   normte = const*normro^(-gamma+1)
normte = proton/(2*kboltz)*normro
if keyword_set(conv) then begin
   rh = rh_s / normro
   T = te_s * normte
   if ~keyword_set(ne_s) then n_e = rh/proton/1.e6 else n_e = ne_s   ; in cgs
   pe = n_e * kboltz * T * 1.e7 ; in cgs
endif else begin
   ; check for CGS
   rh = rh_s
   T = te_s
   if ~keyword_set(ne_s) then n_e = rh/proton else n_e = ne_s
endelse

logT = alog10(T>1.)
logne = alog10(n_e)

sizes=size(rh)
dims=sizes[0]
; if dims<2 exit
if dims eq 3 then begin
   nx=sizes[1]
   ny=sizes[2]
   nz=sizes[3]
endif else begin
   nx = sizes[1]
   nz = sizes[2]
endelse

;if (~(keyword_set(ion))) then begin
;   ion='fe_9'  
;   ion = 'fe_12'
;endif

if (~(keyword_set(nwave))) then begin
   nwave=100
endif
;if ion eq 'fe_9' then w0 = 171.073 ; wave center
;if ion eq 'fe_12' then w0 = 193.509

minwave = w0-0.07
maxwave = w0+0.07
wave=findgen(nwave)/(nwave-1)*(maxwave-minwave)+minwave

pe = n_e * kboltz * T * 1.e7 ; in cgs

ne_sort = sort(n_e)
n_e_sorted = n_e[ne_sort]
Tlg_sorted = logT[ne_sort]

; Read tabulated G(ne,T) values for given number density (n_e_lg) and
; temperature (t_lg) arrays
lookup_goft, ion=ion, w0=w0, dirgot=dirgot,n_e_lg=n_e_lg, logt=t_lg, goft_mat=goft_mat,watom=watom
elements,w0=w0, ion=ion, logTm=logTm, enum=enum, inum=inum, ind=ind

goft=n_e*0.
read_abund,abund_name,abund,abund_ref
line_abunds = abund[enum-1]

; calculate intensity = G(T,n)*ne^2
intens=n_e*0.
for i=0.,n_elements(n_e)-1 do begin
   lc_ne = ([min(abs(n_e_lg-alog10(n_e_sorted[i]))),!c])[1]
   lc_te = ([min(abs(t_lg-tlg_sorted[i])),!c])[1]
   if tlg_sorted[i] lt min(t_lg) or tlg_sorted[i] gt max(t_lg) then goft[ne_sort[i]] = 0. else goft[ne_sort[i]] = interpol(goft_mat[lc_ne,*],t_lg,tlg_sorted[i])
   intens[ne_sort[i]] = goft[ne_sort[i]]*n_e_sorted[i]^2*line_abunds
   print,string(13b)+' % finished: ',float(i)*100./(n_elements(n_e)-1),format='(a,f4.0,$)'
endfor

end
