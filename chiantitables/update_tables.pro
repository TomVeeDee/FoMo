pro update_tables,path=path,spectral_lines=spectral_lines,aia=aia

if (~keyword_set(path)) then path='.'

; find list of all files, and look for aia files (need special treatment) and normal chianti files
listoffiles=file_search(path+'/'+'goft_table_*.dat')
aiafiles=listoffiles[where(strmatch(listoffiles,'*aia*'))]
specfiles=listoffiles[where(~strmatch(listoffiles,'*aia*'))]

if keyword_set(spectral_lines) then begin
;for i=0,n_elements(specfiles)-1 do begin
	for i=0,0 do begin
		; initialise variables for reading
		ion=''
		w0=0.
		print,'Updating ',specfiles[i]
		; for each file, read the ion and rest wavelength
		openr,unit,specfiles[i],/get_lun
		readf,unit,ion
		readf,unit,w0
		free_lun,unit

		; even though all spectroscopic files should be done with sun_coronal.abund (because it is
		; normalised in FoMo is other abundances are needed), let's check if the file name has abco. 
		if (strmatch(specfiles[i],'*abco*')) then goft_table,ion=ion,w0=w0,gotdir=path+'/'
	endfor
endif

if keyword_set(aia) then begin
	make_aiaresponse, sngfilter='all',gotdir=path+'/',file_abund='coronal'
	make_aiaresponse, sngfilter='all',gotdir=path+'/',file_abund='photospheric'
endif

end
