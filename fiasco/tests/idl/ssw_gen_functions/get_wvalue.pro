function get_wvalue, ids
;   
;+
;   Name: get_wvalue 
;
;   Purpose: return widget value 
;
;   Input Parameter:
;      ids - scaler or vector of widget ids
;            (if vector, all values must be ident str)
;
on_error,2
widget_control, ids(0), get_value=retval
sizeval=size(retval)

for i=1,n_elements(ids)-1 do begin
   widget_control, ids(i), get_value=newval
   sizenewval = size(newval)
   if n_elements(sizeval) ne n_elements(sizenewval) then $
      message,'Widget (u)values have incompatible attributes'
   retval=[retval, newval]
endfor
return,retval
end 

