function readgoftcube_dat,emissionsave,data,compress=compress

; reads the data in emissionsave files, not the triangulation
; invoke as data=readgoftcube('fileyouwanttoread')
; data will be a 4xn array, with n the number of grid points in x,y,\lambda

openr,lun,emissionsave,/get_lun,compress=compress
ng=0L
dim=0L
nvars=0L
chiantisize=0LL
readu,lun,dim
readu,lun,ng
readu,lun,nvars
readu,lun,chiantisize
chiantibuffer=bytarr(chiantisize)
readu,lun,chiantibuffer
chiantifile=string(chiantibuffer)
abundsize=0LL
readu,lun,abundsize
abundbuffer=bytarr(abundsize)
readu,lun,abundbuffer
abundfile=string(abundbuffer)
data=fltarr(dim+nvars,ng)
temp=0.
for i=0,dim+nvars-1 do begin
	for j=0,ng-1 do begin
		readu,lun,temp
		data[i,j]=temp
	endfor
endfor

close,lun
free_lun,lun
return,data

end
