;+
; Project     :	SDAC
;
; Name        :	ADD_TAG
;
; Purpose     :	add a tag to a structure
;
; Use         : NEW_STRUCT=ADD_TAG(STRUCT,TAG,TAG_NAME)
;
; Inputs      :	STRUCT = input structure (array or scalar)
;             : TAG_VALUE = tag variable to add
;             : TAG_NAME = tag name
;
; Outputs     :	NEW_STRUCT = new structure
;
; Keywords    :	NAME = new name for structure
;               INDEX = index or tag name where to append new tag [def = last]
;               ERR   = error message [blank if ok]
;               DUPLICATE = set to allow duplicate tag names
;		TOP_LEVEL = If set, then only the top level is searched to
;			    determine if the tag already exists.
;               NO_COPY = do not make copy of input TAG variable
;                         (it will be destroyed after input)
;               NO_PAIR = if adding an array to a structure array, then
;                         add_tag will pair each array element with the
;                         corresponding structure array element. Setting
;                         /NO_PAIR disables this behavior and will add the
;                         entire array to each structure array element.
;
; Restrictions:	Cannot add more than one tag at a time
;
; Category    :	Structure handling
;
; Written     :	Dominic Zarro (ARC)
;
; Version     :	Version 1, 7 November 1994 -- written
;               Version 2, 16 June 1996    -- cleaned up
;		Version 3, 11-Aug-1997, William Thompson, GSFC
;			Added keyword TOP_LEVEL
;		Version 4, 8-Oct-1998, Zarro (SMA/GSFC) - converted to using
;                       CREATE_STRUCT
;               Version 5, 24-Dec-2004, Zarro (L-3Com/GSFC) - vectorized
;                          01-Feb-05, Csillaghy (Univ. Applied Sciences NW Switzerland)
;                          - changed n_elements( struct ) to size(struct /dim ), see at the end.
;                          02-Nov-2005, Kim Tolbert - Previously only handled 1-D structure in
;                            the pairing mode.  Now if structure is > 1-D, and last
;                            dimension matches last dimension of tag_value, then add_tag will
;                            place each tag_value in the corresponding structure element, unless
;                            no_pair is set.
;
;-
;----------------------------------------------------------------------------

function add_tag_err,val,err
err='invalid input'
pr_syntax,'new_struct=add_tag(struct,tag_value,tag_name)'
if exist(val) then return,val else return,-1
end

;----------------------------------------------------------------------------

function add_tag,struct,tag_value,tag_name,index=index,err=err,$
                 duplicate=duplicate,top_level=top_level,no_copy=no_copy,$
                 quiet=quiet,_extra=extra,name=name,no_pair=no_pair

verbose=1-keyword_set(quiet)

;-- bypass for old IDL versions

forward_function add_tag2

if (1-since_version('5.4')) then begin
 if have_proc('add_tag2') then $
  return,add_tag2(struct,tag_value,tag_name,index=index,err=err,$
                 duplicate=duplicate,top_level=top_level,no_copy=no_copy,$
                 quiet=quiet,_extra=extra,name=name)
  message,'no longer supported for this version of IDL - '+!version.release,/cont
  if exist(struct) then return,struct else return,-1
endif

;-- catch errors

err=''
error=0
catch,error
if error ne 0 then begin
 err=err_state()
 if verbose then message,err,/cont
 catch,/cancel
 if exist(struct) then return,struct else return,-1
endif

if (1-exist(tag_value)) then return,add_tag_err(struct,err)

;-- just tag value and name entered

no_copy=keyword_set(no_copy)
new_struct=-1
if  (n_params() eq 2) then begin
 if is_string(tag_value) and exist(struct) then begin
  if no_copy then $
   return,create_struct(tag_value,temporary(struct),name=name) else $
    return,create_struct(tag_value,struct,name=name)
 endif
 return,add_tag_err(struct,err)
endif

;-- no tag name was entered

if is_string(tag_name) then tname=tag_name else tname=''
tname=strtrim(tname(0),2)
if tname eq '' then begin
 err='tag name must be non-blank string'
 if verbose then message,err,/cont
 if exist(struct) then return,struct else return,-1
endif
tname=strupcase(tname)

if n_elements(tname) ne 1 then begin
 err='restricted to adding one tag at a time'
 if verbose then message,err,/cont
 if exist(struct) then return,struct else return,-1
endif

;-- input structure undefined

if (1-is_struct(struct)) and exist(tag_value) then begin
 if no_copy then $
  return,create_struct(tname,temporary(tag_value),name=name) else $
   return,create_struct(tname,tag_value,name=name)
endif

;-- does tag already exist

idl5=idl_release(lower=5)
if (1-keyword_set(duplicate)) or idl5 then begin
 if is_struct(struct) then begin
  if tag_exist(struct,tname,top_level=top_level) then begin
   if idl5 and verbose then message,'duplicate tag - '+tname+' - not added',/cont
   return,struct
  endif
 endif
endif

;-- determine location of added tag

tags=tag_names(struct)
ntags=n_elements(tags)
aindex=ntags-1
if exist(index) then begin
 if is_string(index) then begin
  ilook=where(strup(index) eq tags,icount)
  if icount gt 0 then aindex=ilook[0]
 endif else aindex=index[0]
endif
aindex= -1 > aindex < (ntags-1)
no_copy=keyword_set(no_copy)

;-- check if tag_value is an array with outer dimension
;   equal to array size of input structure. If so, we add each tag_value
;   to each structure by index.

pair=1-keyword_set(no_pair)
ndim_struct = size(struct, /n_dim)

ssiz = size(struct, /dim)
nstruct = last_item(ssiz)

sz=size(tag_value)
nval=sz[n_elements(sz)-3]	;gets size of last dimension
pair=keyword_set(pair)
by_index=0b
if (nstruct gt 1) and (nval eq nstruct) and pair then begin
 by_index=1b
 dim=size(tag_value,/dim)
 if n_elements(dim) ge 2 then begin
  dim=dim[0:n_elements(dim)-2]
 endif else dim=1
 if n_elements(dim) eq 1 then dim=dim[0]
 tvalue=replicate(tag_value[0],dim)
 if n_elements(tvalue) eq 1 then tvalue=tvalue[0]
endif else begin
 if no_copy then tvalue=temporary(tag_value) else tvalue=tag_value
endelse

if aindex eq -1 then new_struct=create_struct(tname,temporary(tvalue))
for i=0,ntags-1 do begin
 if is_struct(new_struct) then $
  new_struct=create_struct(new_struct,tags[i],struct[0].(i)) else $
   new_struct=create_struct(tags[i],struct[0].(i))
 if i eq aindex then new_struct=create_struct(new_struct,tname,temporary(tvalue))
endfor

if is_string(name) then new_struct=create_struct(new_struct,name=strup(name))

orig_dim=size(struct, /dim)
if total(orig_dim) gt 1 then $
 new_struct=replicate2(temporary(new_struct),orig_dim)

struct_assign,struct,new_struct,/nozero

if by_index then begin
 if ndim_struct eq 1 then begin
   if no_copy then new_struct.(aindex+1)=temporary(tag_value) else $
     new_struct.(aindex+1)=tag_value
 endif else begin
   ; we have an array whose last dimension matches the last dimension of the
   ; structure array, but the structure array is more than 1-D.
   ; We have to make tag_value have the same dimensions of the final structure.tag
   ; so we can move it into the structure instead of using a loop.
   ; e.g. if tag_value is (2,4) and struct is (3,4), we want to make
   ; tval (2,3,4) and shove it into structure.tag.  First make
   ; tval (2,1,4) then rebin that to (2,3,4) and it will reproduce the (2,4) 3 times.
   strdim = size(new_struct,/dim)
   ndim = n_elements(strdim)
   newdim = [dim, intarr(ndim-1)+1, strdim[ndim-1] ]
   tval = rebin (reform (tag_value, newdim), size(new_struct.(aindex+1),/dim) )
   if no_copy then new_struct.(aindex+1)=temporary(tval) else $
     new_struct.(aindex+1)=tval
 endelse
endif


return,new_struct & end
