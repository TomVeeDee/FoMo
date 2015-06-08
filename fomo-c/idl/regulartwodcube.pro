pro regulartwodcube,data,xvec,yvec,emiss

; data is the return product of readgoftcube.pro
; it is a 2D array, with 3 columns and N rows
; N = nx*ny
; It is converted into two vectors xvec and yvec (coordinates along the CCD)
; and emiss is the intensity in the CCD (with dimension nx,ny)

; This routines is typically used for imaging data.

sizes=size(data)
N=sizes[2]
temp=where(data[0,*] eq data[0,0],count)
nx=N/count

temp=where(data[1,*] eq data[1,0],count)
ny=N/count

xvec=dblarr(nx)
yvec=dblarr(ny)
emiss=dblarr(nx,ny)

j=0
k=0
l=0
nlambda=1

for i=0,N-1 do begin
	j=i/nlambda
	k=i/nlambda/nx
	j-=k*nx
	xvec[j]=data[0,i]
	yvec[k]=data[1,i]
	emiss[j,k]=data[2,i]
endfor

end
