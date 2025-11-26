FUNCTION mu_file_lines, filename  
;+
; NAME:
;      MU_FILE_LINES
; PURPOSE:
;      Obtain the number of lines in an ascii file
; EXPLANATION:
;      Duplicate of the intrinsic IDL function GET_LINES
;
; CALLING SEQUENCE:
;       MU_GET_LINES,filename
; INPUTS:
;       FILENAME = ASCII file to process
; OUTPUTS:
;       Number of lines in the file
;
; EXAMPLE:
;       NONE
; COMMON BLOCKS:
;       None
; ROUTINES USED:
;
; NOTES
;       None
; MODIFICATION HISTORY:
;       T. Belloni  19 Nov 2008 from http://idlastro.gsfc.nasa.gov/idl_html_help/FILE_LINES.html

   OPENR, unit, filename, /GET_LUN  
   str = ''  
   count = 0ll  
   WHILE ~ EOF(unit) DO BEGIN  
      READF, unit, str  
      count = count + 1  
   ENDWHILE  
   FREE_LUN, unit  
   RETURN, count  
END 
