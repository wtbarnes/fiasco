;+
; Project     : SOHO - CDS
;
; Name        : HAVE_TAG
;
; Category    : Utility
;
; Purpose     :	Checks if structure has a specified tag name
;
; Explanation :	Same as CHKTAG or TAG_EXIST, but looks for
;               substring match
;
; Syntax      : IDL> have=have_tag(struct,tag_name)
;
; Inputs      : STRUCT = structure to check
;               TAG_NAME = tag name to check
;
; Opt. Inputs : None
;
; Outputs     : HAVE = 1 for at least one match, 0 none.
;
; Opt. Outputs: INDEX = index position of found tag
;
; Keywords    : SORT_INDEX = sort found tags indicies
;               EXACT= set for exact match
;               COUNT = # of matches found
;               NO_RECURSE = do not recurse
;               START = tag name match at start of string
;               TAGS = tag names corresponding to index
;               RECURSE = set to recurse
;
; History     : Version 1,  23-Aug-1998,  D.M. Zarro.  Written
;               Version 2,  16-June-2002, Zarro (LAC/GSFC) - added /START
;               Version 3,  22-Dec-2002, Zarro (EER/GSFC) - made NO_RECURSE the default
;
; Contact     : DZARRO@SOLAR.STANFORD.EDU
;-

function have_tag,struct,tag_name,index,sort_index=sort_index,exact=exact,$
                  count=count,no_recurse=no_recurse,start=start,tags=tags,recurse=recurse

have=0b
count=0
sz=size(struct)
stype=sz[n_elements(sz)-2]
sz=size(tag_name)
dtype=sz[n_elements(sz)-2]
index=-1
if (stype ne 8) or (dtype ne 7) then return,have

;-- check for matches

exact=keyword_set(exact)
start=keyword_set(start)
if keyword_set(recurse) then tags=strupcase(rtag_names(struct)) else $
 tags=strupcase(tag_names(struct))
ntags=n_elements(tag_name)

for i=0,ntags-1 do begin
 temp=trim(tag_name[i])
 out=grep(temp,tags,exact=exact,index=tindex,start=start)
 if tindex[0] gt -1 then nindex=append_arr(nindex,tindex)
endfor

;-- sort unique matches

if exist(nindex) then begin
 have=1b
 if keyword_set(sort_index) then nindex=get_uniq(nindex)
 tags=tags(nindex)
endif else begin 
 nindex=-1
 tags=''
endelse

count=n_elements(nindex)
if count eq 1 then begin
 nindex=nindex[0]
 tags=tags[0]
endif
index=nindex

return,have & end
