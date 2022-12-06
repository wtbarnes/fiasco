FUNCTION chk_dir, dir_name, output, fullname=fullname,err=err
;+
; PROJECT:
;       SOHO - CDS
;
; NAME:	
;       CHK_DIR()
;
; PURPOSE:
;       Check the validity of a directory name.
;
; EXPLANATION:
;
; CALLING SEQUENCE: 
;       Result = CHK_DIR(dir_name)
;
; INPUTS:
;       DIR_NAME -- A string specifying the directory name. For VMS system,
;                   a valid directory name can be a logical name, or
;                   any string with a format of '[...]', '[...]aaa.dir', 
;                   or 'aaa.dir'
;
; OPTIONAL INPUTS: 
;       None.
;
; OUTPUTS:
;       RESULT -- 1 if the directory name is valid, 0 otherwise
;
; OPTIONAL OUTPUTS:
;       OUTPUT -- A string indicating the real directory path
;
; KEYWORD PARAMETERS: 
;       FULLNAME -- if set and OUTPUT is present, OUTPUT will contain the full
;                   path specification for the directory
;
; CALLS:
;       DATATYPE, CHKLOG, STR_INDEX, OS_FAMILY
;
; COMMON BLOCKS:
;       None.
;
; RESTRICTIONS: 
;       None.
;
; SIDE EFFECTS:
;       None.
;
; CATEGORY:
;       
; PREVIOUS HISTORY:
;       Written October 9, 1994, by Liyun Wang, GSFC/ARC
;
; MODIFICATION HISTORY:
;       Version 2, Liyun Wang, GSFC/ARC, December 16, 1994
;          Made work for VMS directories
;       Version 3, Liyun Wang, GSFC/ARC, December 29, 1994
;          Added keyword FULLNAME
;          Fixed bug of false reporting if the given string represents 
;             a plain file under Unix OS
;	Version 4, William Thompson, GSFC, 29 August 1995
;		Modified to use OS_FAMILY()
;	Version 5, Zarro, 29 April 1997
;		Added check for blank input
;
; VERSION:
;       Version 5
;-
;
   ON_ERROR, 2
   err=''
   IF N_ELEMENTS(dir_name) EQ 0 THEN BEGIN
      err='Syntax: Result = CHK_PATH(dir_name)'
      return,0
   ENDIF

   IF datatype(dir_name) NE 'STR' THEN BEGIN
      err='Input parameter must be of string type.'
      return,0
   ENDIF

   IF TRIM(dir_name) EQ '' THEN RETURN,0

   output = ''
   IF (!version.os NE 'vms') THEN BEGIN
      IF OS_FAMILY() EQ 'Windows' THEN suffix = '\..' ELSE suffix = '/..'
      IF STRPOS(dir_name,'~') EQ 0 THEN temp = chklog(dir_name) $
      ELSE temp = dir_name
      rr = findfile(temp+suffix)
      IF rr(0) EQ '' THEN RETURN, 0 ELSE BEGIN
         output = temp
         IF KEYWORD_SET(fullname) THEN BEGIN
;---------------------------------------------------------------------------
;           On unix or dos, FINDFILE does not return the full path
;           specification of the directory, we have to do this with two CDs.
;---------------------------------------------------------------------------
            cd, temp, current = old_dir
            cd, old_dir, current = output
         ENDIF
         RETURN, 1
      ENDELSE
      RETURN, 0
   ENDIF
   name_dir = STRUPCASE(dir_name)
   p_dir = STRPOS(name_dir,'.DIR') 
   IF p_dir EQ -1 THEN BEGIN
;----------------------------------------------------------------------
;     DIR_NAME must be either a logical symbol or a directory name
;     with "[..]" format, if it is indeed a directory name
;----------------------------------------------------------------------
      aa = chklog(name_dir)
      IF (aa NE '') THEN BEGIN
         output = aa
         RETURN, 1
      ENDIF
;----------------------------------------------------------------------
;     Try to see if it is already a valid directory name
;----------------------------------------------------------------------
      d_len = STRLEN(name_dir)
      rquote = STRPOS(name_dir,']') 
;----------------------------------------------------------------------
;     Since dir_name is not a logical name, it should have [..] format
;     to be a directory name
;----------------------------------------------------------------------
      IF rquote EQ -1 THEN RETURN, 0

      lquote = STRPOS(name_dir,'[')
      IF lquote EQ -1 THEN BEGIN
         MESSAGE, 'Invalid directory name',/cont
         RETURN, 0
      ENDIF

      IF rquote EQ d_len-1 THEN BEGIN
;----------------------------------------------------------------------
;        directory name in [...] format        
;----------------------------------------------------------------------
         idx = str_index(name_dir,'.')
         IF idx(0) EQ -1 THEN BEGIN
            IF lquote EQ 0 THEN temp_dir = '[000000]'+$
               STRMID(name_dir,1,rquote-1)+'.DIR' $
            ELSE temp_dir = STRMID(name_dir,0,lquote+1)+'000000]'+$
               STRMID(name_dir,lquote+1,rquote-lquote-1)+'.DIR'
         ENDIF ELSE BEGIN
            n_idx = N_ELEMENTS(idx)
            last_name = STRMID(name_dir,idx(n_idx-1)+1,$
                               d_len-idx(n_idx-1)-2)+'.DIR'
            temp_dir = STRMID(name_dir,0,idx(n_idx-1))+']'+last_name
         ENDELSE
      ENDIF ELSE BEGIN
;----------------------------------------------------------------------
;        directory name in [...]xxx format        
;----------------------------------------------------------------------
         temp_dir = name_dir+'.DIR'
      ENDELSE
   ENDIF ELSE temp_dir = dir_name
;----------------------------------------------------------------------
;  DIR_NAME ends with ".dir" and could be a subdirectory
;----------------------------------------------------------------------
   temp = findfile(temp_dir,count = aa)
   f_name = STRCOMPRESS(temp(0))
   IF aa GT 0 THEN BEGIN
      t_len = STRLEN(f_name)
      rquote = STRPOS(f_name,']')
      lquote = STRPOS(f_name,'[')
      root = STRPOS(f_name,'000000')
      p_dir  = STRPOS(f_name,'.DIR')
      IF root EQ -1 THEN $
         output = STRMID(f_name,0,rquote)+'.'+$
         STRMID(f_name,rquote+1,p_dir-rquote-1)+']' $
      ELSE $
         output = STRMID(f_name,0,lquote+1)+$
         STRMID(f_name,root+7,p_dir-rquote-1)+']'
      RETURN, 1
   ENDIF
   RETURN, 0
END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End of 'chk_dir.pro'.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
