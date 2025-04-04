PRO break_file, file, disk_log, dir, filnam, ext, fversion, node, last_dot=last_dot
;+
; Project     : SOHO - CDS
;
; Name        :
;	BREAK_FILE
; Purpose     :
;	Break a filename into its component parts.
; Explanation :
;	Given a file name, break the filename into the parts
;	of disk/logical, the directory, the filename, the
;	extension, and the file version (for VMS)
; Use         :
;	BREAK_FILE, FILE, DISK_LOG, DIR, FILNAM, EXT, FVERSION, NODE
; Inputs      :
;	file	- The file name
; Opt. Inputs :
;	None.
; Outputs     :
;	disk_log- The disk or logical (looks for a ":")
;		  This is generally only valid on VMS machines
;	dir	- The directory
;	filnam	- The filename (excluding the ".")
;	ext	- The filename extension (including the ".")
;	fversion- The file version (only VMS)
;	node	- The Node name (only VMS)
; Opt. Outputs:
;	None.
; Keywords    :
;	LAST_DOT- For filenames which have "." in the name, it takes
;	 	  the LAST "." as the start of the extension.  The
;		  default is from the first ".".
; Calls       :
;	OS_FAMILY
; Common      :
;	None.
; Restrictions:
;	VMS:
;		Assumes that : always precedes []
;	ULTRIX:
;		Right now it has trouble with the ultrix option of use
;		of "." or ".."
; Side effects:
;	None.
; Category    :
;	Utilities, Operating_system.
; Prev. Hist. :
;	Written 1988 by M.Morrison
;	   Aug-91 (MDM) Changed to handle Unix filename convensions
;	28-Feb-92 (MDM) * Adjusted to handle arrays
;	11-Mar-92 (MDM) - Perform a STRTRIM(x,2) on input string before
;			  doing the "break-up"
;	 1-Dec-92 (MDM) - Moved code to do filename, extension and version
;			  number for both VMS and Unix (previously it
;			  did not do version number code for Unix)
;	29-Jan-93 (DMZ/MDM) - checked for node in file name
; Written     :
;	M. Morrison, August 1991.
; Modified    :
;	Version 1, William Thompson, GSFC, 23 April 1993.
;		Incorporated into CDS library.
;	Version 1.1, William Thompson, GSFC, 7 May 1993.
;		Added IDL for Windows compatibility.
;	Version 2, William Thompson, GSFC, 15 June 1995
;		Merged with Yohkoh version.  Added change 11-Nov-93 by D. Zarro
;       	to check for ".]["  and "[000000" in VMS concealed directory
;		names
;	Version 3, William Thompson, GSFC, 29 August 1995
;		Modified to use OS_FAMILY
;       Version 4, Liyun Wang, NASA/GSFC, October 11, 1996
;               Fixed a bug that occurs in MS Windows
;       Version 5, SVH Haugan, UiO, 10 July 1997
;               Made "ifil" loop variable a LONG
;	Version 6, M.Morrison, LMSAL, 5-Nov-97
;		Added /LAST_DOT option
;	Version 7, richard.schwartz@gsfc.nasa.gov, 23-nov-2002
;		protect against '/' in Windows
; Version     :
;	Version 5, 10 July 1997
;
;-
;
   qvms = 1

;--------------------------------------------------------------------------
;  Check (to see if / or \ is in file) if it is VMS (lots of negatives there)
;--------------------------------------------------------------------------
   dummy = WHERE(STRPOS(file, '/') NE -1, count)
   IF (count EQ 0) THEN BEGIN
      dummy = WHERE(STRPOS(file, '\') NE -1, count)
      IF (count NE 0) THEN qvms = 0
   ENDIF ELSE qvms = 0

   n = N_ELEMENTS(file)
   node = STRARR(n)
   disk_log = STRARR(n)
   dir = STRARR(n)
   filnam = STRARR(n)
   ext = STRARR(n)
   fversion = STRARR(n)

   FOR ifil=0L, n-1 DO BEGIN
      file0 = file(ifil)
      file0 = STRTRIM(file0, 2) ;MDM added 11-Mar-92
      len = STRLEN(file0)

;---------------------------------------------------------------------------
;     node name present    ;DMZ added Jan'93
;     (If so then strip it off now and then add it back later)
;---------------------------------------------------------------------------
      dcolon = STRPOS(file0, '::')
      IF dcolon GT -1 THEN BEGIN
         node(ifil) = STRMID(file0, 0, dcolon+2)
         file0 = STRMID(file0, dcolon+2, 1000)
      ENDIF

      IF (qvms) THEN BEGIN
         p = STRPOS(file0, ':')
         IF (p NE 1) THEN disk_log(ifil) = STRMID(file0, 0, p+1) ;includes :
         len = len-p+1
         file0 = STRMID(file0, p+1, len)

;---------------------------------------------------------------------------
;        check for .][ in dir    ;-- DMZ added Nov'93
;---------------------------------------------------------------------------
         IF STRPOS(file0, '.][') NE -1 THEN $
            file0 = str_replace(file0, '.][', '.')

         p = STRPOS(file0, ']')
         IF (p NE -1) THEN dir(ifil) = STRMID(file0, 0, p+1) ;includes ]

;---------------------------------------------------------------------------
;        check for .000000 in dir  ;-- DMZ added Nov'93
;---------------------------------------------------------------------------
         temp = dir(ifil)
         IF STRPOS(temp, '.000000') NE -1 THEN $
            dir(ifil) = str_replace(temp, '.000000', '')
         len = len-p+1
         file0 = STRMID(file0, p+1, len)
      END ELSE IF OS_FAMILY() EQ 'Windows' THEN BEGIN
;---------------------------------------------------------------------------
;        William Thompson, added support for Microsoft Windows, 7 May 1993.
;---------------------------------------------------------------------------
         p = STRPOS(file0, ':')
         IF p NE -1 THEN BEGIN
            disk_log(ifil) = STRMID(file0, 0, p+1) ;Includes :
            len = len - p + 1
            file0 = STRMID(file0, p+1, len)
         ENDIF
         p = -1

;---------------------------------------------------------------------------
;        find last \
;---------------------------------------------------------------------------
         WHILE (STRPOS(file0, '\', p+1) NE -1) DO p = STRPOS(file0, '\', p+1)
         ;Maybe there is a last '/'
         WHILE (STRPOS(file0, '/', p+1) NE -1) DO P = STRPOS(file0, '/', p+1)
         dir(ifil) = STRMID(file0, 0, p+1)
         file0 = STRMID(file0, p+1, len-(p+1))
      END ELSE BEGIN
         p = -1                 ;WTT changed 7-May-93

;---------------------------------------------------------------------------
;       find last /
;---------------------------------------------------------------------------
         WHILE (STRPOS(file0, '/', p+1) NE -1) DO p = STRPOS(file0, '/', p+1)
         dir(ifil) = STRMID(file0, 0, p+1)
         file0 = STRMID(file0, p+1, len-(p+1))
      END

      p = STRPOS(file0, '.')
      if (keyword_set(last_dot)) then p = str_lastpos(file0, '.')
      IF (p EQ -1) THEN BEGIN
         filnam(ifil) = STRMID(file0, 0, len)
         p = len
      END ELSE filnam(ifil) = STRMID(file0, 0, p) ;not include .
      len = len-p
      file0 = STRMID(file0, p, len)
                                ;
      p = STRPOS(file0, ';')
      IF (p EQ -1) THEN BEGIN
         ext(ifil) = STRMID(file0, 0, len)
         p = len
      END ELSE ext(ifil) = STRMID(file0, 0, p) ;includes . but not ;
      len = len-p
      file0 = STRMID(file0, p, len)
                                ;
      fversion(ifil) = ''
      IF (len NE 0) THEN fversion(ifil) = file0

;---------------------------------------------------------------------------
;     now prefix disk name with node name
;---------------------------------------------------------------------------
      IF node(ifil) NE '' THEN disk_log(ifil) = node(ifil)+disk_log(ifil)
   END

   IF (n EQ 1) THEN BEGIN       ;turn output into scalars
      disk_log = disk_log(0)
      dir = dir(0)
      filnam = filnam(0)
      ext = ext(0)
      fversion = fversion(0)
      node = node(0)
   END

END
