pro gaussfitgoftcube,emiss,lambda,int,doppler,sigma,chisq

; this assumes a regular emiss file, where the last dimension coincides with the length of lambda

sizes=size(emiss)
nx=sizes[1]
ny=sizes[2]
nl=n_elements(lambda)

int=dblarr(nx,ny)
doppler=dblarr(nx,ny)
sigma=dblarr(nx,ny)
chisq=dblarr(nx,ny)

l0=mean(lambda)
mem=max(emiss)
errors=mean(emiss)

for i=0,nx-1 do begin
	for j=0,ny-1 do begin
		localmem=max(emiss[i,j,*])
		if (localmem ge mem/100000.) then begin
			yfit=gaussfit(lambda,reform(emiss[i,j,*]),coeff,nterms=3,chisq=ch,estimates=[localmem,l0,(max(lambda)-min(lambda))/4.],measure_errors=replicate(errors,nl))
			int[i,j]=coeff[0]
			doppler[i,j]=(coeff[1]-l0)/l0*2.99e5 ; in km/s
			sigma[i,j]=abs(coeff[2])
			chisq[i,j]=ch
		endif 
	endfor
endfor

end

