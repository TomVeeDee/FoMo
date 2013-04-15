pro regulargoftcube,data,xvec,yvec,lvec,emiss

; data is the return product of readgoftcube.pro
; it is a 2D array, with 4 columns and N rows
; N = nx*ny*nlambda

sizes=size(data)
N=sizes[2]
temp=where(data[0,*] eq data[0,0],count)
nx=N/count

temp=where(data[1,*] eq data[1,0],count)
ny=N/count

temp=where(data[2,*] eq data[2,0],count)
nlambda=N/count

xvec=dblarr(nx)
yvec=dblarr(ny)
lvec=dblarr(nlambda)
emiss=dblarr(nx,ny,nlambda)

j=0
k=0
l=0

for i=0,N-1 do begin
	j=i/nlambda
	k=i/nlambda/nx
	j-=k*nx
	l=i mod nlambda 
	;if ((j mod nlambda) eq 0) then xvec[j]=data[0,i]
	;if ((k mod nlambda*nx) eq 0) then yvec[k]=data[1,i]
	;if ((j eq 0) and (k eq 0)) then lvec[l]=data[2,i]
	xvec[j]=data[0,i]
	yvec[k]=data[1,i]
	lvec[l]=data[2,i]
	emiss[j,k,l]=data[3,i]
endfor

end
