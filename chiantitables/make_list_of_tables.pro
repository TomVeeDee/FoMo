pro make_list_of_tables,path=path

if (~keyword_set(path)) then path='.'

; find list of all files, and look for aia files (need special treatment) and normal chianti files
listoffiles=file_search(path+'/'+'goft_table_*.dat')
aiafiles=listoffiles[where(strmatch(listoffiles,'*aia*'))]
specfiles=listoffiles[where(~strmatch(listoffiles,'*aia*') and ~strmatch(listoffiles,'*dkist*'))]

openw,out,'list_of_tables.txt',/get_lun
printf,out,'# element',string(9B),'wavelength (A)',string(9B),'minlevel',string(9B),'maxlevel',string(9B),'filename'

for i=0,n_elements(specfiles)-1 do begin
	ion=''
	w0=0.
	atweight=0
	minlevel=0
	maxlevel=0
	; for each file, read the ion and rest wavelength
	openr,unit,specfiles[i],/get_lun
	readf,unit,ion
	readf,unit,w0
	readf,unit,atweight
	readf,unit,minlevel,maxlevel
	free_lun,unit

	printf,out,ion,string(9B),w0,string(9B),minlevel,string(9B),maxlevel,string(9B),specfiles[i]
endfor
free_lun,out

end
