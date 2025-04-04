function event_name,event
;
;+
;-
event_name=tag_names(event,/structure_name)
return, strmid(event_name,7,strlen(event_name))
end
