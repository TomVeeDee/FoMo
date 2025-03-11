function readgoftcube_txt,emissionsave,data,compress=compress

; reads the data in emissionsave files, not the triangulation
; invoke as data=readgoftcube('fileyouwanttoread')
; data will be a 4xn array, with n the number of grid points in x,y,\lambda

openr,lun,emissionsave,/get_lun,compress=compress
header=lonarr(3)
readf,lun,header
stringheader=strarr(2)
readf,lun,stringheader
data=fltarr(header[0]+header[2],header[1])
readf,lun,data

close,lun
free_lun,lun
return,data

end
