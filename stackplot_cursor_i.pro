pro STACKPLOT_CURSOR_I, stackplot=stackplot, min_time=min_time, max_time=max_time, $
 			line_x=line_x, line_y=line_y, vel=vel, d_vel=d_vel, plot=plot, no_transpose=no_transpose, quiet=quiet, soho=soho, color=color, linestyle=linestyle, ds=ds, dt=dt, text_out=text_out, text_color=text_color, boxcar=boxcar, line_width=line_width, line_pix=line_pix, line_int = line_int, line_d_int = line_d_int

;PURPOSE
;  To measure velocities from plot_st_stackplot
;  by point & click using mouse cursor
;
;INPUT
;  None
;
;OPTIONAL INPUT
;Transpose			- (keyword) has to be set if the plot_st_stackplot used it.
;				  Othwerise returns 1/v instead of v
;SOHO				- (keyword) set if the stackplots are generated from SOHO data
;				  or other instrument at L1; since then the arc second is slightly shorter
;Quiet				- (keyword) set if you don't want any console output
;
;OUTPUT
;  None
;
;OPTIONAL OUTPUT
;  Line_x, line_y		- Line coordinates (for plot_st_stackplot_line, TBD)
;  Vel				- calculated velocities
;  D_vel			- velocity uncertainty (crude estimate)
;  DS, DT			- uncertainties in the spatial and temporal resolution (default 0.75 arcsec and 6 s)
;  STACKPLOT			- variable containing the stackplot if intensities are to be extracted
;  BOXCAR			- smoothing in the STACKPLOT_CURSOR
;  LINE_PIX			- how many pixels in the calculation of intensities
;
;PROGRAMMING NOTES
;  1) Note that the program expects to be used after STACKPLOT_CURSOR.
;  2) Note that the keyword REVERSE does not need to be set; as the cursor command is used with /data.
;
;  3) If intensities along the line are to be outputted, the program requires setting STACKPLOT (along with MIN_TIME, MAX_TIME, BOXCAR, etc.) 
;     the same way as STACKPLOT_CURSOR
;  
;
;HISTORY
;2015-12-17	JD		- adapted from map_line_cursor.pro
;2016-01-11	JD		- added calculation of the velocity uncertainty based on a 'rod estimate'
;				  with rod width being given by AIA resolution (1.5")
;				  (longest path in shortest time vs shortest path in logest time in a 'rod' on stackplot)
;2016-01-12	JD		- formatted the output for line_x and line_y to be copy-paste ready
;				- added DS, DT as keywords
;2016-01-20	JD		- changed calculation of D_vel and D_vel0 using the sqrt sigma formula
;2018-01-08	JD		- deprecated TRANSPOSE; added NO_TRANSPOSE; corresponding to changes in STACKPLOT_CURSOR
;				- implemented REVERSE
;2018-03-13	JD		- added keyword TEXT_OUT, TEXT_COLOR
;2018-03-27	JD		- modified text output (commented out some lines)
;2022-04-21	JD		- changed default color from 0 to 255
;2025-03-20	JD		- corrected a bug in velocity calculation with /NO_TRANSPOSE
;
;2025-06-19	JD		- started adding calculation of intensities along the line
;				- implemented BOXCAR (set if plot_st_stackplot was used)
;				- LINE_PIX - if the line should have finite width
;
;2025-07-01	JD		- added MIN_TIME and MAX_TIME
;				- disabled the -1 in the  NTI = n_elements(tind);-1
;				  (why was it even there?!)
;2025-07-09	JD		- work on intensity extraction (unfinished)
;
;2025-07-10	JD		- CHANGED DS and DT to mean DS and DT in stackplot, not uncertainties of the click ! (which are now ds/2 and dt/4)
;				- first prototype of intensity extraction finished
;				  cleaned up a bit
;				- UNTESTED line_pix averaging !
;2025-07-17	JD		- added back-conversion of min_time, max_time at the end of the program
;
;2025-08-12	JD		- changed calculation of grid points for /LINE_PIX



if 	keyword_set(soho)	then arcsec2km = 718. else arcsec2km = 725.
if not(keyword_set(color))	then color=255
if not(keyword_set(linestyle))	then linestyle=2




; I.
; Housekeeping - stackplot time coordinate manipulation
; (this taken over en bloc from PLOT_ST_STACKPLOT)
;------------------------------------------------------
if keyword_set(stackplot) then begin

  
  scoord	= stackplot.scoord
  tcoord	= stackplot.tcoord
  int_st	= stackplot.stackplot
  
  if not(keyword_set(DS)) then ds = stackplot.ds
  if not(keyword_set(DT)) then dt = stackplot.dt
  ; ds		= stackplot.ds
  ; dt		= stackplot.dt

  ;min_time & max_time OK?
  if keyword_set(min_time) then begin
    if not(is_string(min_time)) then begin
      print, '-!-  STACKPLOT_CURSOR:  MIN_TIME must be a string'
      stop
    endif

    min_time = anytim(min_time, out_style='vms')
    if (strmid(min_time, 0, 11) EQ ' 1-Jan-1979') then strput, min_time, strmid(tcoord[0], 0, 11), 0
    if (min_time LT min(tcoord)) or (min_time GT max(tcoord)) then begin
 	print, (min_time LT min(tcoord)), '  ', min_time, '  ', min(tcoord)
 	print, (min_time GT max(tcoord)), '  ', min_time, '  ', max(tcoord)
	print, '-!-  STACKPLOT_CURSOR:  INCORRECT time range: MIN_TIME'
	print, ''
; 	stop
    endif
  endif
  if keyword_set(max_time) then begin
    if not(is_string(max_time)) then begin
      print, '-!-  STACKPLOT_CURSOR:  MAX_TIME must be a string'
      stop
    endif

    max_time = anytim(max_time, out_style='vms')
    if (strmid(max_time, 0, 11) EQ ' 1-Jan-1979') then strput, max_time, strmid(tcoord[n_elements(tcoord)-1], 0, 11), 0
    if (max_time GT max(tcoord)) or (max_time LT min(tcoord)) then begin
  	print, (max_time GT max(tcoord)), '  ', max_time, '  ', max(tcoord)
	print, (max_time LT min(tcoord)), '  ', max_time, '  ', min(tcoord)
	print, '-!-  STACKPLOT_CURSOR:  INCORRECT time range: MAX_TIME'
	print, ''
; 	stop
    endif
  endif

  if keyword_set(trange) and (keyword_set(min_time) or keyword_set(max_time)) then begin
	print, '-!-  STACKPLOT_CURSOR:  Use either trange or min_time/max_time keywords'
	stop
  endif
  if keyword_set(trange) then begin
	min_time	= min(trange) ;-dt/2d
	max_time	= max(trange) ;+dt/2d
  endif

  if keyword_set(interval) and keyword_set(ytickinterval) then begin
	print, '-!-  STACKPLOT_CURSOR:  Both INTERVAL and YTICKINTERVAL are set. Using INTERVAL only'
	ytickinterval	= interval
  endif
  if keyword_set(interval) and not(keyword_set(ytickinterval)) then ytickinterval = interval


  ; if not(keyword_set(ds))		then ds		= (deriv(scoord))[0]
  if not(keyword_set(srange)) 		then srange	= [min(scoord)-ds/2d, max(scoord)+ds/2d]
  if not(keyword_set(min_time)) 	then min_time	= time_px_shift(min(tcoord), dt, /subtract)
  if not(keyword_set(max_time)) 	then max_time	= time_px_shift(max(tcoord), dt, /add)
  ; if not(keyword_set(trange)) 	then trange	= [ max([min(tcoord), min_time]), min([max(tcoord), max_time]) ]
  ; if not(keyword_set(trange)) 	then trange	= [ max([min(tcoord)-dt/2d, min_time]), min([max(tcoord)+dt/2d, max_time]) ]
 

  ;we no longer need trange, min_time and max_time should now be properly defined.
  ;note that time offsets in the tcoord are probably not properly handled
  ;e.g., does the time correspond to the center of the 'pixel' as for the spatial coordinate?

  ;truncate the scoord and tcoord to include only the points within the specified ranges
  sind		= where((scoord -ds/2d GE srange[0]) and (scoord +ds/2d LE srange[1]))
  tind		= where((time_px_shift(tcoord, dt, /subtract) GE min_time) $
		    and (time_px_shift(tcoord, dt, /add)      LE max_time))
  ; tind		= where((tcoord -dt/2d GE min_time)  and (tcoord +dt/2d LE max_time))
  scoord	= scoord[sind]
  tcoord	= tcoord[tind]
  NXI		= n_elements(sind)
  NTI		= n_elements(tind);-1
  
  ; help, int_st, sind, tind, int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1]]

  if keyword_set(true) and (!D.name EQ 'PS') then begin 
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
  endif
  if not(keyword_set(true)) or (!D.name EQ 'X') then begin
    dimension_array	= [1,0]
    if keyword_set(boxcar) then begin
	int_st	= smooth(int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1]], boxcar, nan=nan, edge_truncate=edge_truncate)
    endif else begin
	int_st	= int_st[sind[0]:sind[NXI-1],tind[0]:tind[NTI-1]]
    endelse
  endif

  ; if keyword_set(reverse) then begin
  ; ;   srange = reverse(srange)
  ; ;   scoord	= reverse(scoord)
  ;   if keyword_set(no_transpose) then begin
  ; ;     int_st = reverse(int_st)
  ;     int_st = reverse(int_st, 1)
  ;     xrange = [srange[1], srange[0]]
  ;   endif else begin
  ; ;     int_st = transpose(reverse(int_st, 2), dimension_array)
  ;     int_st = reverse(int_st, 1)
  ;     yrange = [srange[1], srange[0]]
  ;   endelse
  ; endif else begin
  ;   if keyword_set(no_transpose) then begin
  ;      xrange = srange
  ;   endif else begin
  ;      yrange = srange
  ;   endelse
  ; endelse

  ttcoord	= anytim2utc(tcoord)
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
	print, '-!-  STACKPLOT_CURSOR:  Too many days!', Ndays
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
  ; time		= ttcoord.time/1d3 
  time		= ttcoord.time/1d3  -tmin
  
  
  NS		= NXI
  NT		= NTI
  
endif 





; II.
; Here we produce & plot the line
; (clicking on the plotted stackplot)
;------------------------------------
jumpBegin:
delvarx, line_x, line_y

!mouse.button = 0
i=0	& x=dblarr(100) 	& y=dblarr(100)

while (!mouse.button NE 4) do begin
;   print, !mouse.button
  if (!mouse.button EQ 2) then begin
    goto, jumpBegin 
  endif else begin
    cursor, xc, yc, /data
    x[i]=xc & y[i] = yc
    plots, x[i], y[i], psym=1, /data
    if (i GE 1) then plots, x[i-1:i], y[i-1:i], /data, color=color, linestyle=linestyle
    i=i+1 & wait, 0.2
  endelse
endwhile
line_X=x[0:i-2]
line_Y=y[0:i-2]

; if keyword_set(vertical) then line_X[1:-1] = line_X[0]

if     keyword_set(plot)	then plots, line_X, line_Y, /data
if not(keyword_set(quiet))	then begin
  out_line_x	= '; line_x = ['
  out_line_y	= '; line_y = ['
  for j=0, i-3, 1 do begin
    out_line_x	= out_line_x +string(line_x[j], format='(F9.2)')+', '
    out_line_y	= out_line_y +string(line_y[j], format='(F9.2)')+', '
  endfor
  out_line_x	= out_line_x +string(line_x[i-2], format='(F9.2)')+']'
  out_line_y	= out_line_y +string(line_y[i-2], format='(F9.2)')+']'
  print, out_line_x
  print, out_line_y
  print, ''
endif


vel		= dblarr(i-2)			; velocity
D_vel		= dblarr(i-2)			; first-order uncertainty
D_vel0		= dblarr(i-2)			; first-order uncertainty based on ds=0.75


if not(keyword_set(ds)) then ds = 1.5		; [arcsec, AIA]
if not(keyword_set(dt)) then dt = 12.		; [s, AIA]

;calculate velocities
;Note: x=scoord, y=tcoord. If /transpose is set, we have to do ^(-1)
for j=0, i-3, 1 do begin

  ; vel[j]	= (line_x[j+1] -line_x[j]) / (line_y[j+1] -line_y[j]) / arcsec2km				; km s^(-1)
  ; if not(keyword_set(no_transpose)) then vel[j] = 1d/vel[j]

  if keyword_set(no_transpose) then begin
    vel[j] = (line_x[j+1] -line_x[j]) * arcsec2km / (line_y[j+1] -line_y[j]) 
  endif else begin
    vel[j] = (line_y[j+1] -line_y[j]) * arcsec2km / (line_x[j+1] -line_x[j]) 
  endelse 
    
    
  ;uncertainty. Take the position uncertainty to be 1/4 the AIA resolution, i.e., 1.5"/4
  ;             Take the time uncertainty to be 1/4 the cadence, i.e., 3 s
  ; as in Lorincik et al. 2025
  ;
  ; this changed here as of 2025-07-10
  
  if    (keyword_set(no_transpose)) then begin
;     D_vel[j]	= sqrt(  arcsec2km^2d *2d* ds^2d *(1d/(line_y[j+1]-line_y[j]))^2d 			$
; 			+2d*dt^2d *((line_x[j+1]-line_x[j])*arcsec2km/(line_y[j+1]-line_y[j])^2d)^2d )
;     D_vel0[j]	= sqrt(  arcsec2km^2d *2d* (0.75d)^2d *(1d/(line_y[j+1]-line_y[j]))^2d 			$
; 			+2d*6d^2d *((line_x[j+1]-line_x[j])*arcsec2km/(line_y[j+1]-line_y[j])^2d)^2d )
    D_vel[j]	= sqrt(  arcsec2km^2d *2d* (ds/4.)^2d *(1d/(line_y[j+1]-line_y[j]))^2d 			$	
			+2d*(dt/4.)^2d *((line_x[j+1]-line_x[j])*arcsec2km/(line_y[j+1]-line_y[j])^2d)^2d )
    D_vel0[j]	= sqrt(  arcsec2km^2d *2d* (1.5/4.)^2d *(1d/(line_y[j+1]-line_y[j]))^2d 			$
			+2d*3.^2d *((line_x[j+1]-line_x[j])*arcsec2km/(line_y[j+1]-line_y[j])^2d)^2d )
 endif else begin
;     D_vel[j]	= sqrt(  arcsec2km^2d *2d* ds^2d *(1d/(line_x[j+1]-line_x[j]))^2d 			$
; 			+2d*dt^2d *((line_y[j+1]-line_y[j])*arcsec2km/(line_x[j+1]-line_x[j])^2d)^2d )
;     D_vel0[j]	= sqrt(  arcsec2km^2d *2d* (0.75d)^2d *(1d/(line_x[j+1]-line_x[j]))^2d 			$
; 			+2d*6d^2d *((line_y[j+1]-line_y[j])*arcsec2km/(line_x[j+1]-line_x[j])^2d)^2d )
    D_vel[j]	= sqrt(  arcsec2km^2d *2d* (ds/4.)^2d *(1d/(line_x[j+1]-line_x[j]))^2d 			$
			+2d*(dt/4.)^2d *((line_y[j+1]-line_y[j])*arcsec2km/(line_x[j+1]-line_x[j])^2d)^2d )
    D_vel0[j]	= sqrt(  arcsec2km^2d *2d* (1.5/4.)^2d *(1d/(line_x[j+1]-line_x[j]))^2d 			$
			+2d*3.^2d *((line_y[j+1]-line_y[j])*arcsec2km/(line_x[j+1]-line_x[j])^2d)^2d )
  endelse

  if not(keyword_set(quiet)) 	then print, '; Velocity  for segment ', trim(j+1), ' is ',		$
					       string(vel[j], format='(F9.2)'), ' km/s   +- ',		$
					       string(D_vel[j], format='(F9.2)'), ' km/s    (', 	$
					       string(D_vel0[j], format='(F9.2)'), ' km/s)'

					       
					              
  ; DISABLED FOR NOW - UNFINISHED AS OF 2025-06-19, 2025-06-30, 2025-07-01
  ; NOTES TBD:
  ; - tcoord conversion for stackplot interpolation  		-- DONE
  ; - xap, yap as function of j (for different line segments)	-- DONE, UNTESTED
  ; - intensity average and uncertainty				-- DONE
  ; - test with real stackplot					-- DONE

  ; intensities along line segments
  ; (this portion of code is used/modified from st_stackplot.pro)
   
  if keyword_set(stackplot) then begin  
  
   
    sgn	 	= +1d
    if (line_X[j] NE line_X[j+1]) then begin

	  ;non-vertical cuts
	  a		= (line_Y[j+1] -line_Y[j]) / (line_X[j+1] -line_X[j])
	  b		=  line_Y[j] -a*line_X[j]
	  dst		= sqrt( (line_X[j] -line_X[j+1])^2d +(line_Y[j] -line_Y[j+1])^2d )
	  ; DS		= 1d /sqrt( ((cos(atan(a)))/DX)^2d +((sin(atan(a)))/DY)^2d ) 	; this is already known/given
	  ; npx		= fix(dst /DS +0.5d) +1
	  
          ; DS NE DT, which means we will have to vary step size depending on the line slope:
	  if (atan(abs(a)) GT atan(DS/DT)) then begin
	    delta 	= DS /sin(atan(abs(a)))
	    omega	= DS /cos(atan(abs(a)))
    	  endif else begin
    	    delta	= DT /cos(atan(abs(a)))
    	    omega	= DS /cos(atan(abs(a)))
          endelse
	  
	  npx		= fix(dst/delta +0.5) +1
	  phi		= replicate(atan(a), npx)
	  

	  ; ; horizontal case
   ;        ; if (line_Y[j] EQ line_Y[j+1]) and (n_elements(line_x) EQ 2) then begin
   ;        if (line_Y[j] EQ line_Y[j+1]) then begin
	  ;   lcoord 	= dindgen(npx) *DS +line_X[j]
	  ; endif else begin
	  ;   lcoord	= dindgen(npx) *DS
	  ; endelse

	  if (line_x[j+1] LT line_x[j]) then sgn=-1d
	  xp		= sgn *delta *dindgen(npx) *cos(phi[j]) +line_X[j]			; if tan(phi[j]) = a, phi = <0,!dpi/2)
	  yp		= sgn *delta *dindgen(npx) *sin(phi[j]) +line_Y[j]
    endif else begin

	  ;vertical cuts
	  ; (this does not make sense for stackplots, but leave it here consistency)
	  npx		= fix(abs(Line_Y[j+1] -Line_Y[j]) /DS +0.5d) +1
	  omega		= dT
	  
	  if (line_y[j+1] LT line_Y[j]) then sgn = -1d
	  xp		= replicate(line_X[j], npx)
	  yp		= sgn *DS *dindgen(npx) +line_Y[j]
	  
	  ; if keyword_set(straight_cut) or (n_elements(line_y) EQ 2) then begin	  
	  ;   lcoord	= DS *dindgen(npx) +line_Y[j]
	  ; endif else begin
	  ;   lcoord	= DS *dindgen(npx)
	  ; endelse

  	  phi		= replicate(sgn *!dpi/2d, npx)
    endelse

    ; Is there finite a slit width?
    ; if yes, then define the slit (aperture)
    ; and calculate everything in pixels across the cut
    if keyword_set(line_width) then line_pix = fix(line_width/delta +0.5)

    if keyword_set(line_pix) then begin

  	line_width	= line_pix *DS
  	print, '% STACKPLOT_CURSOR: Line_width corrected to '+trim(line_width)+' arc sec'

  	xap		= dblarr(npx,2*line_pix+1)
  	yap		= dblarr(npx,2*line_pix+1)
  
      for ipx=0, npx-1, 1 do begin
       for jpx=-line_pix, +line_pix, 1 do begin
	xap[ipx, jpx+line_pix]	= xp[ipx] +jpx*omega*sin(phi[i])
	yap[ipx, jpx+line_pix]	= yp[ipx] -jpx*omega*cos(phi[i])					; minus?
        endfor
      endfor
    endif else begin

  	xap		= xp
  	yap		= yp
    endelse

    ; for output
    if NOT(keyword_set(line_pix)) then line_width = 0d

    NSLIT_PIX	= n_elements(xap[0,*])
    int_line	= dblarr(npx,NSLIT_PIX)

    for jpx=0, NSLIT_PIX-1, 1 do begin
  	; if not(keyword_set(boxcar)) then begin 

  	  int_line[*,jpx]= interpolate(int_st, $
				       interpol(dindgen(NS), scoord, yap[*,jpx]), $
				       interpol(dindgen(NT), time,   xap[*,jpx]))
  	; endif else begin
          ; 
  	  ; int_line[*,jpx]= interpolate(smooth(stackplot.stackplot, boxcar, nan=nan, edge_truncate=edge_truncate), $
				 ;      interpol(dindgen(NT), tcoord, 	      xap[*,jpx]), $
				 ;      interpol(dindgen(NS), stackplot.scoord, yap[*,jpx]) )   
  	; endelse
    endfor
  
    if not(keyword_set(quiet))	then begin
    
 	print, '; Intensity for segment '+trim(j+1) +' is '+string(average(int_line), format='(F9.2)') +' DN' + $
    	       '   +- '  +string(stddev(int_line), format='(F9.2)') +' DN'
    
    	; now average over the slit width if necessary
    	if (NSLIT_PIX GT 1) then begin
    	  out_int_line 	= total(int_line, 2, /double)/double(NSLIT_PIX)
    	  print, '; Intensity for segment '+trim(j+1) +' is '+string(average(out_int_line), format='(F9.2)') +' DN' + $
    	         '   +-  ' +string(stddev(out_int_line), format='(F9.2)')+' DN'
        endif else begin
          out_int_line 	= average(int_line)
          out_int_line_e = stddev(int_line)
        endelse

      line_int = out_int_line
      line_d_int = out_int_line_e

    	; str_int_line	= ';  int_line = ['
    	; for ipx=0, npx-2, 1 do begin
    	;   str_int_line= str_int_line +trim(string(out_int_line[ipx], format='(F9.2)'))+', '
    	; endfor
    	; str_int_line	= str_int_line +trim(string(out_int_line[ipx], format='(F9.2)'))+']'
  
    	; print, str_int_line
      ;   print, '; '
    endif
  endif  
  

endfor




; III.
; plot the line and text (velocities)
if keyword_set(text_out) then begin

  if n_elements(line_x) EQ 2 then begin
    if not(is_string(text_out)) then begin
      if (fix(text_out) EQ 1) then txt	= trim(string(vel[0], format='(F9.1)'))+'+-'+trim(string(D_vel[0], format='(F9.1)'))+' km/s'
    endif else begin
      txt = text_out 
    endelse
    
    a 		= linfit(line_x[0:1]/(!x.crange[1]-!x.crange[0])*!d.x_size, line_y[0:1]/(!y.crange[1]-!y.crange[0])*!d.y_size)
    orientation	= atan(a[1])*180d/!dpi
    sign	= (-1d)*(orientation LT 0d) + (+1d)*(orientation GE 0d)
    xyouts, line_x[0] +(!x.crange[1]-!x.crange[0])/1d3, line_y[0]+(!y.crange[1]-!y.crange[0])/1d2, txt, /data, $
	    orientation=orientation, color=text_color

  endif
endif



; convert back
if keyword_set(max_time) then begin
  min_time	= anytim(min_time, out_style='vms')
  max_time	= anytim(max_time, out_style='vms')
endif  

end
