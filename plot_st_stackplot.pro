function time_px_shift, time, dt, add=add, subtract=subtract
; shifts the time by dt/2
  tt	= (anytim2utc(time))
  if keyword_set(add)      then tt.time=tt.time+dt*1d3/2d
  if keyword_set(subtract) then tt.time=tt.time-dt*1d3/2d
  out	= anytim(tt, out_style='vms')
return, out
end




;====================================================================================
PRO plot_st_stackplot, st_stackplot, srange=srange, trange=trange, min_time=min_time, max_time=max_time, dmin=dmin, dmax=dmax, log=log, channel=channel, true=true, isotropic=isotropic, charsize=charsize, charthick=charthick, thick=thick, xstyle=xstyle, ystyle=ystyle, xmargin=xmargin, ymargin=ymargin, xticks=xticks, yticks=yticks, xminor=xminor, yminor=yminor, xtickname=xtickname, ytickname=ytickname, xtickinterval=xtickinterval, ytickinterval=ytickinterval, interval=interval, xticklen=xticklen, yticklen=yticklen, xtitle=xtitle, ytitle=ytitle, title=title, inverse=inverse, no_erase=no_erase, erase=erase, no_transpose=no_transpose, reverse=reverse, seconds=seconds, boxcar=boxcar, nan=nan, edge_truncate=edge_truncate, bgcs=bgcs, bgts=bgts, consm=consm, xedge=xedge, yedge=yedge, sob=sob, sharp=sharp, cmask=cmask    ;, t_offset=t_offset

;INPUT
;--------------------
;st_stackplot			- stackplot structure (output from st_stackplot.pro)
;				  tags: stackplot - (2D array) [int units]	stackplot
;					ds	  - (double)   [arcsec]		spatial resolution
;					dt	  - (double)   [time, VMS]	temporal resolution
;					scoord	  - (1D array) [arcsec]		spatial coordinates (center of pixels)
;					tcoord    - (1D array) [time, VMS]	temporal coordinates (center of pixels)
;

;PROGRAMMING NOTES
;--------------------
;scoord				- always corresponds to the CENTER of the spatial pixel
;tcoord				- TB IMPLEMENTED
;
;
;
;MODIFICATION HISTORY
;--------------------
;2014-08-05	JD		- added keyword DS and implemented half-pixel shifts in Scoord
;				- implemented /erase and /inverse
;				- implemented /transpose
;
;2014-11-21	JD		- changed input to st_stackplot only.
; 				- added min_time and max_time; will influence trange
;				  rewritten the handling of indices throughout the program
;
;				- added the TIME_PX_SHIFT function
;2014-12-01	JD		- changed /erase to /no_erase: /erase is kept for compatibility, but not functional
;				- added conversion of min_time and max_time to VMS style.
;				- added KEYWORDS xtickinterval and ytickinterval
;				  added handling of time ticks if axis parameters not set.
;				- corrected missing 1d3 factor in TIME_PX_SHIFT
;
;2014-12-10	JD		- added keyword SECONDS: prevents automatic ticks in time format, shows 'actual' time in seconds
;				  Useful if measurements of speeds etc. are to be made.
;2015-05-13	JD		- Changed time axis label if  /seconds is set
;2015-12-11	JD		- Changed the default ytickinterval for better plotting of a longer time series
;2015-12-17	JD		- Minor changes to error messages
;2016-01-12	JD		- added keyword BOXCAR for smoothing the stackplot
;				- implemented REVERSE to reverse the spatial direction
;2016-01-14	JD		- disabled stop of the program if min_time or max_time incorrect
;2016-10-31	JD		- added keyword CHANNEL (for true-color plotting)
;				- corrected a bug involving midnight in the data
;2016-11-17	JD		- corrected handling of YTICKNAME (was useless because of the
;				  if (yticks EQ 0) and NOT(keyword_set(seconds)) then condition,
;				  after which the yticknames were defined
;2016-11-29	JD		- added keyword INTERVAL (same meaning as YTICKINTERVAL but less confusing)
;2017-04-05	JD		- added "NOT(keyword_set(yticks))" to 'automatic tick interval calculation'
;2017-04-26	JD		- added keyword TRUE - for PS device, similar implementation as for plotmap.pro
;				  note this required adding a variable "dimension_array" for TRANSPOSE
;2017-09-18	JD		- added keywords XSTYLE, YSTYLE
;				  If not set, the program re-sets them to 1 (as before with /xs, /ys in plot)
;				- added keywords XTICKLEN, YTICKLEN
;
;2017-09-18	JD		CHANGED:
;				- TRANSPOSE to NO_TRANSPOSE
;				- [hopefully] changed all respective keywords (i.e., xminor, etc. now refers to X axis
; 									 independently of what is on it)
;
;2017-09-26	JD		CORRECTED:
;				- plotting with REVERSE
;2018-01-25	JD		- added /NAN keyword: if BOXCAR is set and there are NaN or Inf values, use /NAN
;				- added /EDGE_TRUNCATE to SMOOTH call
;
;2018-02-05	JD		- changed automatic setting for DMIN and DMAX (dmin=1 if /log); taken from plotmap.pro
;				- added additional checks for MIN_TIME and MAX_TIME
;				  (... or min_time GT max(tcoord); max_time LT min(tcoord) )
;				- added handling of MIN_TIME and MAX_TIME with times only (does not require date)
;2018-02-06	JD		- added checks for min_time, max_time being strings
;2018-02-22	JD		- Original MIN_TIME and MAX_TIME are now stored and returned (so that the input keyword is not overwritten)

;2019-11-25     JL              - added a bunch of keywords which can be used for post-processing of the input stackplot
;                               - USAGE:

;                               - First, use image filters based on convolution, esentially controlled by two keywords: consm and cmask
;                               - consm: before the image is convolved with an arbitrary kernel, it is smoothed with a boxcar
;                               - cmask: when set, the images from the convolution and the original ones are normalized, added, and normalized again
;                               - sharp: kernel for making the image sharper, has very similar effect to the Laplacian function
;                               - xedge: kernel for detection of edges along the Y direction 
;                               - yedge: kernel for detection of edges along the X direction 
;                               - sobel: combines the latter two using a built-in function in IDL

;                               - One can then use some background subtraction: 
;                               - bgts: for each position along the cut, the BG averaged along all times is subtracted 
;                               - bgcs: for each time, the BG averaged along the cut is subtracted 
    

scoord		= st_stackplot.scoord
tcoord		= st_stackplot.tcoord
int_st		= st_stackplot.stackplot
ds		= st_stackplot.ds
dt		= st_stackplot.dt
; dt		= 0d


;min_time & max_time OK?
if keyword_set(min_time) then begin
  if not(is_string(min_time)) then begin
    print, '-!-  PLOT_ST_STACKPLOT:  MIN_TIME must be a string'
    stop
  endif

  min_time = anytim(min_time, out_style='vms')
  if (strmid(min_time, 0, 11) EQ ' 1-Jan-1979') then strput, min_time, strmid(tcoord[0], 0, 11), 0
  if (min_time LT min(tcoord)) or (min_time GT max(tcoord)) then begin
; 	print, ''
	print, '-!-  PLOT_ST_STACKPLOT:  INCORRECT time range: MIN_TIME'
	print, ''
; 	stop
  endif
endif
if keyword_set(max_time) then begin
  if not(is_string(max_time)) then begin
    print, '-!-  PLOT_ST_STACKPLOT:  MAX_TIME must be a string'
    stop
  endif

  max_time = anytim(max_time, out_style='vms')
  if (strmid(max_time, 0, 11) EQ ' 1-Jan-1979') then strput, max_time, strmid(tcoord[n_elements(tcoord)-1], 0, 11), 0
  if (max_time GT max(tcoord)) or (max_time LT min(tcoord)) then begin
; 	print, ''
	print, '-!-  PLOT_ST_STACKPLOT:  INCORRECT time range: MAX_TIME'
	print, ''
; 	stop
  endif
endif

if keyword_set(trange) and (keyword_set(min_time) or keyword_set(max_time)) then begin
	print, '-!-  PLOT_ST_STACKPLOT:  Use either trange or min_time/max_time keywords'
	stop
endif
if keyword_set(trange) then begin
	min_time	= min(trange) ;-dt/2d
	max_time	= max(trange) ;+dt/2d
endif

if keyword_set(interval) and keyword_set(ytickinterval) then begin
	print, '-!-  PLOT_ST_STACKPLOT:  Both INTERVAL and YTICKINTERVAL are set. Using INTERVAL only'
	ytickinterval	= interval
endif
if keyword_set(interval) and not(keyword_set(ytickinterval)) then ytickinterval = interval

; if not(keyword_set(ds))		then ds		= (deriv(scoord))[0]
if not(keyword_set(srange)) 	then srange	= [min(scoord)-ds/2d, max(scoord)+ds/2d]
if not(keyword_set(min_time)) 	then min_time	= time_px_shift(min(tcoord), dt, /subtract)
if not(keyword_set(max_time)) 	then max_time	= time_px_shift(max(tcoord), dt, /add)
; if not(keyword_set(trange)) 	then trange	= [ max([min(tcoord), min_time]), min([max(tcoord), max_time]) ]
; if not(keyword_set(trange)) 	then trange	= [ max([min(tcoord)-dt/2d, min_time]), min([max(tcoord)+dt/2d, max_time]) ]

if not(keyword_set(charthick))	then charthick	= 1 +(!D.name EQ 'PS')
if not(keyword_set(thick))	then thick	= 1 +(!D.name EQ 'PS')
if not(keyword_set(charsize))	then charsize	= 1.25
if not(keyword_set(thick))	then thick	= 2
; if not(keyword_set(dmin))	then dmin	= min(int_st)
; if not(keyword_set(dmax))	then dmax	= max(int_st)
if not(keyword_set(dmin))	then dmin	= max([min(int_st), keyword_set(log)])
if not(keyword_set(dmax))	then dmax	= min([max(int_st), 1d5])
if not(keyword_set(xstyle))	then xstyle	= 1
if not(keyword_set(ystyle))	then ystyle	= 1
if not(keyword_set(xticklen))	then xticklen	= 0.01
if not(keyword_set(yticklen))	then yticklen	= 0.01
if not(keyword_set(xmargin))	then xmargin	= [10,3]
if not(keyword_set(ymargin))	then ymargin	= [4,1.5]
if not(keyword_set(xticks))	then xticks	= 0
if not(keyword_set(yticks))	then yticks	= 0
if not(keyword_set(xminor))	then xminor	= 0
if not(keyword_set(yminor))	then yminor	= 0
if not(keyword_set(xtickname))	then xtickname	= replicate('', xticks+1)
if not(keyword_set(ytickname))	then ytickname	= replicate('', yticks+1)
if not(keyword_set(xtitle))	then xtitle	= 'Position along cut [arc sec]'
if not(keyword_set(ytitle))	then begin
 if keyword_set(seconds) 	then ytitle	= 'Time from the '+min_time+' [s]' else ytitle	= 'Time [UT]'
endif
if not(keyword_set(title))	then title	= ''

if not(keyword_set(consm))	then consm	= 1

;add xtickinterval, ytickinterval here?

if keyword_set(xedge) or keyword_set(yedge) or keyword_set(sharp) then begin 

    proc = int_st
    
    kernel = fltarr(3,3)

    if keyword_set(xedge) then kernel = [[-1, 0, 1],[-2, 0, 2],[-1, 0, 1]]      ;kernels depend on various functions
    if keyword_set(yedge) then kernel = [[ 1, 2, 1],[ 0, 0, 0],[-1,-2,-1]]
    if keyword_set(sharp) then kernel = [[-1,-1,-1],[-1, 8,-1],[-1,-1,-1]]
    
    print, kernel 
    
    proc = convol(smooth(proc, [consm,consm]), kernel)
    proc = shift(proc, 0, 3)                            ; smoothes the shift of the processed image, which probably originates in the 3x3 conv. kernel 
    
    if keyword_set(cmask) then begin 
    
        proc    = proc+abs(min(proc))                   ; make it positive 
        nproc   = proc/max(proc)
        norig   = int_st/max(int_st)
        new_st  = nproc + norig
        new_st  = new_st/max(new_st)                    ; normalization to 1.0
        proc    = new_st
        
        print, 'The output stackplot has been normalized'

    endif 
    
    int_st = proc 
    
endif 

if keyword_set(sob) then begin

    proc=sobel(int_st)
    
    if keyword_set(cmask) then begin 
    
        proc    = proc+abs(min(proc))                   
        nproc   = proc/max(proc)
        norig   = int_st/max(int_st)
        new_st  = nproc + norig
        new_st  = new_st/max(new_st)                 
        proc    = new_st
        
        print, 'The output stackplot has been normalized'

    endif 
    
    int_st = proc 

endif 


if keyword_set(bgts) then begin

    proc = int_st

    nc=n_elements(proc(*,0))    
    means=fltarr(nc)
    
    for i=0, nc-1 do begin           
    
        means[i] = mean(proc(i,*), /nan)          
        proc[i,*] = proc[i,*] - means[i]   
        
    endfor 
    
    int_st = proc
        
endif 

if keyword_set(bgcs) then begin

    proc = int_st
    
    nt=n_elements(proc(0,*))
    means=fltarr(nt)
  
    for i=0, nt-1 do begin                
    
        means[i] = mean(proc(*,i), /nan)          
        proc[*,i] = proc[*,i] - means[i]   
        
    endfor 
    
    int_st = proc
    
endif 


;we no longer need trange, min_time and max_time should now be properly defined.
;note that time offsets in the tcoord are probably not properly handled
;e.g., does the time correspond to the center of the 'pixel' as for the spatial coordinate?

;truncate the scoord and tcoord to include only the points within the specified ranges
sind		= where((scoord -ds/2d GE srange[0]) and (scoord +ds/2d LE srange[1]))
tind		= where((time_px_shift(tcoord, dt, /subtract) GE min_time) $
		    and (time_px_shift(tcoord, dt, /add) LE max_time))
; tind		= where((tcoord -dt/2d GE min_time)  and (tcoord +dt/2d LE max_time))
scoord		= scoord[sind]
tcoord		= tcoord[tind]
NXI		= n_elements(sind)
NTI		= n_elements(tind)-1

if keyword_set(true) then begin 
  dimension_array	= [1,0,2]			; for no_transpose
  if keyword_set(boxcar) then begin
; 	int_st	= smooth(int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1],*], boxcar)
	int_st	= int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1],*]
	for ist=0, 2 do begin
	  int_st[*,*,ist]	= smooth(reform(int_st[*,*,ist]), boxcar, nan=nan, edge_truncate=edge_truncate)
	endfor
  endif else begin
	int_st	= int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1],*]
  endelse
endif else begin
  dimension_array	= [1,0]
  if keyword_set(boxcar) then begin
	int_st	= smooth(int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1]], boxcar, nan=nan, edge_truncate=edge_truncate)
  endif else begin
	int_st	= int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1]]
  endelse
endelse

if keyword_set(reverse) then begin
;   srange = reverse(srange)
;   scoord	= reverse(scoord)
  if keyword_set(no_transpose) then begin
;     int_st = reverse(int_st)
    int_st = reverse(int_st, 1)
    xrange = [srange[1], srange[0]]
  endif else begin
;     int_st = transpose(reverse(int_st, 2), dimension_array)
     int_st = reverse(int_st, 1)
   yrange = [srange[1], srange[0]]
  endelse
endif else begin
  if keyword_set(no_transpose) then begin
     xrange = srange
  endif else begin
     yrange = srange
  endelse
endelse




ttcoord		= anytim2utc(tcoord)
min_time	= anytim2utc(min_time)
max_time	= anytim2utc(max_time)

;correct for midnights in datasets from various days
;untested as of 2014-11-21
Ndays		= (max_time.mjd - min_time.mjd)
if (Ndays GT 0) then begin
	max_time.mjd	= min_time.mjd
	max_time.time	= max_time.time +(3600L *1000L *24L *long(Ndays))
endif
if (Ndays GT 9) then begin
	print, ''
	print, '-!-  PLOT_ST_STACKPLOT:  Too many days!', Ndays
	print, ''
	stop
endif
if (Ndays GT 0) then begin
  for iday=1L, long(Ndays), 1L do begin
    inextday			= where(ttcoord.mjd EQ ttcoord[0].mjd +iday)
    ttcoord[inextday].mjd	= ttcoord[0].mjd
    ttcoord[inextday].time	= ttcoord[inextday].time +(3600L *1000L *24L *long(iday))
  endfor
endif

;set tmin and tmax, in seconds
tmin		= min(min_time.time)/1d3
tmax		= max(max_time.time)/1d3
time		= ttcoord.time/1d3


;set up the plot
if NOT(keyword_set(no_erase)) then erase
if keyword_set(inverse) then begin
  C = -1d
  dmax2 = dmin
  dmin2 = dmax
endif else begin
  C = 1d
  dmax2 = dmax
  dmin2 = dmin
endelse


;automatic time tick intervals
IF keyword_set(no_transpose) THEN BEGIN
  if NOT(keyword_set(yticks)) and NOT(keyword_set(ytickinterval)) then begin
    ytickinterval=60
    if (tmax-tmin) GE 1200d    then ytickinterval=600 
    if (tmax-tmin) GE 7200d    then ytickinterval=1200
    if (tmax-tmin) GE 3600d *6 then ytickinterval=3600
  endif
  if (yticks EQ 0) and NOT(keyword_set(seconds)) and NOT(keyword_set(ytickname)) then begin
    ntimeint	= fix((tmax-tmin)/double(ytickinterval))
    ytickname	= strarr(ntimeint+1)
    for it=0, ntimeint, 1 do begin
      if ((min_time.time +it*double(ytickinterval)*1d3) / (3600L *1000L *24L)) LT 1 then begin
        ytickname[it] 	= anytim({mjd:min_time.mjd,      time:min_time.time +it*double(ytickinterval)*1d3} , out_style='vms')
      endif else begin
        days		= (min_time.time +it*double(ytickinterval)*1d3) / (3600L *1000L *24L)
        ytickname[it]	= anytim({mjd:min_time.mjd+days, time:min_time.time +it*double(ytickinterval)*1d3 - (3600L *1000L *24L *long(days))}, out_style='vms')
      endelse
        ytickname[it] 	= strmid(ytickname[it], 12, 5)
    endfor
  endif
ENDIF ELSE BEGIN
  if NOT(keyword_set(xticks)) and NOT(keyword_set(xtickinterval)) then begin
    xtickinterval=60
    if (tmax-tmin) GE 1200d    then xtickinterval=600 
    if (tmax-tmin) GE 7200d    then xtickinterval=1200
    if (tmax-tmin) GE 3600d *6 then xtickinterval=3600
  endif
  if (xticks EQ 0) and NOT(keyword_set(seconds)) and NOT(keyword_set(xtickname)) then begin
    ntimeint	= fix((tmax-tmin)/double(xtickinterval))
    xtickname	= strarr(ntimeint+1)
    for it=0, ntimeint, 1 do begin
      if ((min_time.time +it*double(xtickinterval)*1d3) / (3600L *1000L *24L)) LT 1 then begin
       xtickname[it] 	= anytim({mjd:min_time.mjd,      time:min_time.time +it*double(xtickinterval)*1d3} , out_style='vms')
      endif else begin
        days		= (min_time.time +it*double(xtickinterval)*1d3) / (3600L *1000L *24L)
        xtickname[it]	= anytim({mjd:min_time.mjd+days, time:min_time.time +it*double(xtickinterval)*1d3 - (3600L *1000L *24L *long(days))}, out_style='vms')
      endelse
        xtickname[it] 	= strmid(xtickname[it], 12, 5)
    endfor
  endif
ENDELSE

;plot
if keyword_set(no_transpose) then begin

	plot, 	/nodata, /norm, [srange[0], srange[1]], [0d, tmax-tmin], 		$
				xrange=xrange,						$
				xstyle=xstyle, ystyle=ystyle,				$
		charsize=charsize, charthick=charthick, /noerase,			$
		thick=thick, xthick=thick, ythick=thick, isotropic=isotropic,		$
		xmargin=[xmargin[0],xmargin[1]], ymargin=[ymargin[0],ymargin[1]], 	$
		xticklen=-xticklen, yticklen=-yticklen, xticks=xticks, yticks=yticks,	$
		xminor=xminor, yminor=yminor, xtickname=xtickname, ytickname=ytickname, $
		xtickinterval=xtickinterval, ytickinterval=ytickinterval,		$
		xtitle=xtitle, ytitle=ytitle, title=title

	if NOT(keyword_set(reverse)) then begin
	  lower_left	= convert_coord(min(scoord) -ds/2d, min(time)-tmin -dt/2d, /data, /to_device)
	  upper_right	= convert_coord(max(scoord) +ds/2d, max(time)-tmin +dt/2d, /data, /to_device)
	endif else begin
	  lower_left	= convert_coord(max(scoord) +ds/2d, min(time)-tmin -dt/2d, /data, /to_device)
	  upper_right	= convert_coord(min(scoord) -ds/2d, max(time)-tmin +dt/2d, /data, /to_device)
	endelse

	if (!D.name EQ 'PS') then begin
	if (keyword_set(log)) then begin
		tv, bytscl(C*(alog10(int_st[*,*,*] > (dmin) < (dmax))),					$
			min=(C*alog10(dmin2)), max=(C*alog10(dmax2)), /nan), true=true,			$
		lower_left(0,0), lower_left(1,0), /device,						$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endif else begin
		tv, bytscl(C*(int_st[*,*,*] > (dmin) < (dmax)),						$
			min=(C*(dmin2)), max=(C*(dmax2)), /nan), true=true,				$
		lower_left(0,0), lower_left(1,0), /device,						$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endelse
	endif
	if (!D.name EQ 'X') then begin
	if (keyword_set(log)) then begin
		tv, bytscl(congrid(C*(alog10(int_st > (dmin) < (dmax))),				$
				upper_right(0,0)-lower_left(0,0)+1, upper_right(1,0)-lower_left(1,0)+1),$
			min=(C*alog10(dmin2)), max=(C*alog10(dmax2)), /nan),  channel=channel,		$
		lower_left(0,0), lower_left(1,0), /device,						$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endif else begin
		tv, bytscl(congrid(C*(int_st > (dmin) < (dmax)),					$
				upper_right(0,0)-lower_left(0,0)+1, upper_right(1,0)-lower_left(1,0)+1),$
			min=(C*(dmin2)), max=(C*(dmax2)), /nan), channel=channel,			$
		lower_left(0,0), lower_left(1,0), /device,						$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endelse
	endif


	plot, 	/nodata, /norm, [srange[0], srange[1]], [0d, tmax-tmin], 		$
				xrange=xrange,						$
				xstyle=xstyle, ystyle=ystyle,				$
		charsize=charsize, charthick=charthick, /noerase,			$
		thick=thick, xthick=thick, ythick=thick, isotropic=isotropic,		$
		xmargin=[xmargin[0],xmargin[1]], ymargin=[ymargin[0],ymargin[1]], 	$
		xticklen=xticklen, yticklen=yticklen, xticks=xticks, yticks=yticks,	$
		xminor=xminor, yminor=yminor,						$
		xtickinterval=xtickinterval, ytickinterval=ytickinterval,		$
		xtickname=replicate(' ', xticks+20), ytickname=replicate(' ', yticks+20)
endif else begin

	plot, 	/nodata, /norm, [0d, tmax-tmin], [srange[0], srange[1]], 		$
				yrange=yrange,						$
				xstyle=xstyle, ystyle=ystyle,				$
		charsize=charsize, charthick=charthick, /noerase,			$
		thick=thick, xthick=thick, ythick=thick, isotropic=isotropic,		$
		xmargin=[xmargin[0],xmargin[1]], ymargin=[ymargin[0],ymargin[1]], 	$
		xticklen=-xticklen, yticklen=-yticklen, xticks=xticks, yticks=yticks,	$
		xminor=xminor, yminor=yminor, xtickname=xtickname, ytickname=ytickname, $
		xtickinterval=xtickinterval, ytickinterval=ytickinterval,		$
		xtitle=ytitle, ytitle=xtitle, title=title

	if NOT(keyword_set(reverse)) then begin
	  lower_left	= convert_coord(min(time)-tmin -dt/2d, min(scoord) -ds/2d, /data, /to_device)
	  upper_right	= convert_coord(max(time)-tmin +dt/2d, max(scoord) +ds/2d, /data, /to_device)
	endif else begin
	  lower_left	= convert_coord(min(time)-tmin -dt/2d, max(scoord) +ds/2d, /data, /to_device)
	  upper_right	= convert_coord(max(time)-tmin +dt/2d, min(scoord) -ds/2d, /data, /to_device)
	endelse

	if (!D.name EQ 'PS') then begin
	if (keyword_set(log)) then begin
		tv, bytscl(C*(alog10(transpose(int_st[*,*,*], dimension_array) > (dmin) < (dmax))),	$
			min=(C*alog10(dmin2)), max=(C*alog10(dmax2)), /nan),  				$
		lower_left(0,0), lower_left(1,0), /device, true=true,					$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endif else begin
		tv, bytscl(C*(transpose(int_st[*,*,*], dimension_array) > (dmin) < (dmax)),		$
			min=(C*(dmin2)), max=(C*(dmax2)), /nan), 					$
		lower_left(0,0), lower_left(1,0), /device, true=true,					$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endelse
	endif
	if (!D.name EQ 'X') then begin
	if (keyword_set(log)) then begin
		tv, bytscl(congrid(C*(alog10(transpose(int_st, dimension_array) > (dmin) < (dmax))),	$
				upper_right(0,0)-lower_left(0,0)+1, upper_right(1,0)-lower_left(1,0)+1),$
			min=(C*alog10(dmin2)), max=(C*alog10(dmax2)), /nan),  channel=channel,		$
		lower_left(0,0), lower_left(1,0), /device,						$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endif else begin
		tv, bytscl(congrid(C*(transpose(int_st, dimension_array) > (dmin) < (dmax)),		$
				upper_right(0,0)-lower_left(0,0)+1, upper_right(1,0)-lower_left(1,0)+1),$
			min=(C*(dmin2)), max=(C*(dmax2)), /nan), channel=channel,			$
		lower_left(0,0), lower_left(1,0), /device,						$
		xsize=upper_right(0,0)-lower_left(0,0)+1, ysize=upper_right(1,0)-lower_left(1,0)+1
	endelse
	endif

	plot, 	/nodata, /norm, [0d, tmax-tmin], [srange[0], srange[1]], 		$
				yrange=yrange,						$
				xstyle=xstyle, ystyle=ystyle,				$
		charsize=charsize, charthick=charthick, /noerase,			$
		thick=thick, xthick=thick, ythick=thick, isotropic=isotropic,		$
		xmargin=[xmargin[0],xmargin[1]], ymargin=[ymargin[0],ymargin[1]], 	$
		xticklen=xticklen, yticklen=yticklen, xticks=xticks, yticks=yticks,	$
		xminor=xminor, yminor=yminor,						$
		xtickinterval=xtickinterval, ytickinterval=ytickinterval,		$
		xtickname=replicate(' ', xticks+20), ytickname=replicate(' ', yticks+20)
endelse


min_time	= anytim(min_time, out_style='vms')
max_time	= anytim(max_time, out_style='vms')
end
