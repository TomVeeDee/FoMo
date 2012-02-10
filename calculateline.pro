pro calculateline, rho, T, wave, ion, goft

; input: 	rho	mass density in kg/m^-3
;		T	temperature in K
; 		wave	wavelength range (Angstrom)
;		ion	selected ion (e.g. 'fe_9')
; output:	goft	G(T,n) in the wavelength range (erg cm^3 sr^-1 s^-1)

ioneq_name= !xuvtop+'/ioneq/chianti.ioneq'
abund_name=!xuvtop+'/abundance/sun_photospheric.abund'

proton=1.67262158*10^(-27.)

n_e=rho/proton/10^6.

if (n_elements(wave) eq 0) then begin
	minwave=160.
	maxwave=180.
	ion='fe_9'
	nwave=100
	wave=findgen(nwave)/(nwave-1)*(maxwave-minwave)+minwave
endif

ch_synthetic,min(wave),max(wave),output=line,density=n_e,logt_isothermal=alog10(T),ioneq_name=ioneq_name,sngl_ion=ion

make_chianti_spec, line, wave, out, instr_fwhm=0.01, abund_name=abund_name

goft=out.spectrum

end
