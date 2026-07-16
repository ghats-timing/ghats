pro mu_fxbread,unit,dato,colonna,riga

common dati,tempo,evento,first_time   ; common block for keeping the data

if(first_time eq 1) then begin
	; If this is the first call, read in the whole table (limited to the
	; relevant columns)
	; notice that data type of evento is automatically recognized
	fxbreadm,unit,[1,2],tempo,evento

; Now the data are available in memory. Give out the desired one
    case colonna of
		1: dato = tempo(riga)
		2: dato = evento(riga)
	endcase
end