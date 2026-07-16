pro gh_help,alpha=alpha,task=task,help=help
;+
; NAME:
;      GH_HELP
; PURPOSE:
;      Open the curated GHATS routine directory.
; EXPLANATION:
;      GH_HELP displays one of the curated HELP/ghats_routines_*.txt
;      files from inside GHATS/IDL through PYTHON/gh_help.
;      By default it opens the task-organized section. Use /ALPHA for the
;      alphabetical section. Inside the pager, search with /text.
;
; CALLING SEQUENCE:
;      GH_HELP
;      GH_HELP,/TASK
;      GH_HELP,/ALPHA
;      GH_HELP,/HELP
;
; KEYWORDS:
;      TASK  = show the task-organized section. This is the default.
;      ALPHA = show the alphabetical section.
;      HELP  = print this usage message.
;
; NOTES:
;      The pager defaults to less. From IDL, GH_HELP opens the pager in a
;      separate terminal window so that less receives normal keyboard input.
;      Set GH_HELP_TERMINAL to choose a terminal emulator, for instance xterm.
;      Set GH_HELP_PAGER to override the pager, for instance GH_HELP_PAGER=print
;      for non-interactive tests. Set GH_HELP_PYTHON to choose a specific
;      Python executable if direct execution is not available. Set GH_HELP_FONT
;      and GH_HELP_FONT_SIZE to override the terminal font. The same
;      functionality is available from the shell with:
;         gh_help
;         gh_help alpha
;         gh_help ,alpha
;-

if(keyword_set(help)) then begin
   print,''
   print,'GH_HELP'
   print,''
   print,'Open the GHATS routine directory.'
   print,''
   print,'Usage:'
   print,'  gh_help'
   print,'  gh_help,/task'
   print,'  gh_help,/alpha'
   print,''
   print,'Inside the pager, search with /text, for example /coherence.'
   print,'Shell: gh_help  or  gh_help alpha  or  gh_help ,alpha'
   print,'Env:   GH_HELP_TERMINAL selects terminal; GH_HELP_PAGER selects pager'
   print,'       GH_HELP_PYTHON selects Python; GH_HELP_FONT selects font'
   print,''
   return
endif

if(keyword_set(alpha) and keyword_set(task)) then begin
   print,'ERROR: Use only one of /ALPHA or /TASK.'
   retall
endif

section = 'task'
if(keyword_set(alpha)) then begin
   section = 'alpha'
endif

root = getenv('MU_PATH')
if(root eq '') then root = '.'
script = root + '/PYTHON/gh_help'

if(~file_test(script)) then begin
   script = 'PYTHON/gh_help'
endif

if(~file_test(script)) then begin
   print,'ERROR: Cannot find PYTHON/gh_help.'
   print,'Set MU_PATH to the GHATS repository root or run from that directory.'
   retall
endif

python = getenv('GH_HELP_PYTHON')
if(python eq '') then begin
   spawn,script + ' ' + section + ' terminal'
endif else begin
   spawn,python + ' ' + script + ' ' + section + ' terminal'
endelse

end
