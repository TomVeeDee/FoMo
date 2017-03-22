function readgoftcube,emissionsave,data

fileparts=strsplit(emissionsave,'.',/extract)

compress=0
if (fileparts[-1] eq 'gz') then begin
	compress=1
endif

case fileparts[-1-compress] of 
	'txt': data=readgoftcube_txt(emissionsave,compress=compress)
	'dat': data=readgoftcube_dat(emissionsave,compress=compress)
	else: print, 'Extension ', fileparts[-1-compress],' not recognised'
endcase

return,data

end
