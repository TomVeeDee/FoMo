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

for i=0,nx-1 do begin
	for j=0,ny-1 do begin
		yfit=gaussfit(lambda,reform(emiss[i,j,*]),coeff,nterms=4,chisq=ch,estimates=[max(emiss[i,j,*]),l0,(max(lambda)-min(lambda))/4.,min(emiss[i,j,*])])
		int[i,j]=coeff[0]
		doppler[i,j]=(coeff[1]-l0)/l0*2.99e5 ; in km/s
		sigma[i,j]=abs(coeff[2])
		chisq[i,j]=ch
	endfor
endfor

end

