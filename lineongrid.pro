pro lineongrid, rho, T, emission,outwave,wave=wave,nwave=nwave,minwave=minwave,maxwave=maxwave,ion=ion

; input: 	rho		array of mass density in kg/m^-3
;		T		array of temperature in K
; output: 	emission	emission as a function of x, y, z (array[nx,ny,nz,nwave]
;		outwave		the wavelength range of the spectroscopic direction in emission

ioneq_name= !xuvtop+'/ioneq/chianti.ioneq'
abund_name=!xuvtop+'/abundance/sun_photospheric.abund'

proton=1.67262158*10^(-27.)

n_e=rho/proton/10^6.

sizes=size(rho)
dims=sizes[0]
; if dims<2 exit
nx=sizes[1]
ny=sizes[2]
nz=sizes[3]

if (~(keyword_set(nwave))) then begin
	nwave=100
endif
if (~(keyword_set(minwave))) then begin
	minwave=170.9
endif
if (~(keyword_set(maxwave))) then begin
	maxwave=171.2
endif

if (~(keyword_set(wave))) then begin
	wave=findgen(nwave)/(nwave-1)*(maxwave-minwave)+minwave
endif
if (~(keyword_set(ion))) then begin
	ion='fe_9'
endif
outwave=wave

if (dims eq 2) then begin
	emission=dblarr(nx,ny,nwave)
	nz=1
endif
if (dims eq 3) then begin
	emission=dblarr(nx,ny,nz,nwave)
endif

; CAREFUL: for a 100x150x150 data cube, this will approximately take 24 hours

; since ch_synthetic sorts the T, we'd better also sort n_e so we have control

t_sort=sort(reform(alog10(T),n_elements(T)))
n_e_sorted=n_e[t_sort]
T_sorted=T[t_sort]

ch_synthetic,min(wave),max(wave),output=allline,density=reform(n_e_sorted,n_elements(n_e)),logt_isothermal=reform(alog10(T_sorted),n_elements(T)),ioneq_name=ioneq_name,sngl_ion=ion,/no_sum_int

; the above command sorts T!!!

ch_synthetic,min(wave),max(wave),output=singleline,density=n_e[0,0,0],logt_isothermal=alog10(T[0,0,0]),ioneq_name=ioneq_name,sngl_ion=ion

; find out how to convert the information in 'line' to a real spectrum, i.e. do something like make_chianti_spec.pro
; perhaps do a for loop over the all the results from ch_synthetic?

for i=0,nx*ny*nz-1 do begin
	singleline.logt_isothermal[0]=allline.logt_isothermal[i]
	singleline.lines.int[0]=allline.lines.int[i]
	make_chianti_spec, singleline, wave, out, abund_name=abund_name
	coords=array_indices(rho,t_sort[i])
	if (dims eq 2) then begin
		emission[coords[0],coords[1],*]=out.spectrum
	endif
	if (dims eq 3) then begin
		emission[coords[0],coords[1],coords[2],*]=out.spectrum
	endif
endfor



end
