pro gaussfitgoftcube,emiss,lambda,peak=peak,doppler=doppler,sigma=sigma,chisq=chisq,int=int

; this assumes a regular emiss file, where the last dimension coincides with the length of lambda (given in \AA)
; peak is the peak intensity of the spectral line (erg cm^-2 s^-1 sr^-1 \AA^-1, if the wavelength is given in \AA)
; int is the integrated intensity of the _fitted_ spectral line (erg cm^-2 s^-1 sr^-1)
; doppler is the fitted Doppler shift (in km/s)
; sigma is the width of the fitted Gaussian (in \AA)
; chisq is the chi-squared of the fit, and indicates the non-Gaussianity of the spectral line.

sizes=size(emiss)
nx=sizes[1]
ny=sizes[2]
nl=n_elements(lambda)

peak=dblarr(nx,ny)
int=dblarr(nx,ny)
doppler=dblarr(nx,ny)
sigma=dblarr(nx,ny)
chisq=dblarr(nx,ny)

l0=mean(lambda)
mem=max(emiss)
errors=mean(emiss)

for i=0,nx-1 do begin
	for j=0,ny-1 do begin
		localmem=max(emiss[i,j,*],maxpos)
		if (localmem ge mem/100000.) then begin
			yfit=gaussfit(lambda,reform(emiss[i,j,*]),coeff,nterms=3,chisq=ch,estimates=[localmem,lambda[maxpos],(max(lambda)-min(lambda))/4.],measure_errors=replicate(errors,nl))
			peak[i,j]=coeff[0] ; coeff[0] is the peak intensity, because it is the constant in front of the Gaussian
			doppler[i,j]=(coeff[1]-l0)/l0*2.99e5 ; in km/s
			sigma[i,j]=abs(coeff[2])
			chisq[i,j]=ch
			; to go from the peak intensity to total intensity
			; we take that y = A\exp(-x^2/2\sigma^2), with A is coeff[0], and \sigma is coeff[2]
			; Then A = I/\sigma/\sqrt{2\pi}
			int[i,j]=coeff[0]*abs(coeff[2])*sqrt(2*!Pi)
		endif 
	endfor
endfor

end

