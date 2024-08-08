;PURPOSE
;create a ST_stackplot: stack images over a slit

;INPUT
;map			- map array
;cut_X			- X position of the CUT [2-element vector if /straight_cut] (optionally named slit_x)
;cut_Y			- Y position of the CUT [2-element vector if /straight_cut] (optionally named slit_y)

;OPTIONAL INPUT
;min_time, max_time	- time constraints. If either not set, the beginning/end of the
; 			  map time series will be assumed.
;straight_cut		- if set, the Cut_X and Cut_Y should be 2-element vectors
;			  if they are not, this program will cut the vector into 2 elements only
;
;slit_width / slit_pix	- half-width of the slit across the cut
;			  -!- ST_STACKPLOT: only one of these keywords should be set
;			  SLIT_WIDTH [arc sec]			- half-width of the slit
;			  SLIT_PIX   [pixels, in DS units]	- the width of the slit is 2*slit_pix +1 in px units
;
;
;boxcar			- smoothes AIA data with the given boxcar value before producing the stackplot
;/SPLINE		- if set, splines over the cut nodes to create a smooth curve


;OUTPUT (RETURN)
;out.stackplot		- resulting stackplot; a 2D array [S,T]
;out.tcoord		- spatial coordinate [date, time]
;out.scoord		- spatial coordinate [arc sec]

;NOTES
;			1)  The user should consider whether to use the /SPLINE or not
;			    If /SPLINE is not set, the program calculates stackplot along a piece-wise linear cut
;
;			2a) The program calculates its own spacing across the cut - variable DS
;			    In principle for maps where (DX NE DY), the DS will depend on the direction of the first cut segment.
;			    The DS is then set constant for the rest of the cut.
;			2b) The location of the nodes 2, ..., n of the cut is modified by this program to conform to the DS value
;			    This modification is less than a half-pixel
;
;			3)  The user should consider the width of the slit (SLIT_WIDTH or SLIT_PIX) carefully
;			    This is true especially for cuts that have acute angles,
;			    where large slits will lead to interpolation over the same image pixels.
;
;

;HISTORY
;2014-03-10	JD 	- written
;2014-08-05	JD	- corrected scoord if slit is in the Solar X direction
;2014-10-19	JD	- added option to input Slit_X and Slit_Y directly
; 			  added keyword /straight_slit: if set, the Slit_X and Slit_Y
; 			    should be 2-element vectors as was expected till now
;			- rewritten the handling of maps spanning different dayap
;2014-11-21	JD	- added tag DT: temporal resolution of the plot
;			  calculated as the median of the time resolution
;2014-12-04	JD 	- coordinates now as a function of the map index;
;			  this has to be done as some maps in a series
;			  may be shifted (e.g., have different XC)
;			- added keyword FAST to skip this.
;2015-12-16	JD	- added XP and YP to the OUT : slit positions
;			  This should be useful for plotting the slit axis on maps
;			- added cut_x, cut_y, but kept slit_x, slit_y
;2016-01-12	JD	- added keyword BOXCAR to smooth the AIA data before producing the stackplot
;			  Note that this keyword is NOT the same as boxcar in plot_st_stackplot,
; 			  which smoothes the final stackplot only.
;2016-10-31	JD	- corrected a bug involving midnight in the data
;2017-09-26	JD	- changed interval=0.6 in spline_p to interval=(DX+DY)/2d
;
;2018-01-18	JD	- REWRITE OF THE STACKPLOT CALCULATION
;
; 			- added XP and YP as output keywords (not allowed as inputs)
;			- added SLIT_WIDTH and SLIT_PIX
;			  changed atan(a) to phi and set phi = !dpi/2d for vertical cuts
;			  this produces a ~1e-17 value of cos(phi[0]) for YS calculation, but this should be OK (?)
;			- changed the way npx is calculated: not fix(dst/DS) +1, but fix(dst/DS +0.5) +1
;			- added keyword SPLINE - if spline_p is to be used on the cut points
;			- changed the condition for reverse: cut_x[n_elements(cut_x)-1] instead of cut_x[1]
;			- the keyword straight_cut now implies that the cut has only 2 nodes, this program discards the rest
;2018-01-22	JD	(contd.)
;			- fixed various bugs
;			- implementation of slit_width and slit_pix in progress
;			- changed output format of various messages to include 'ST_STACKPLOT' program name
;			- added cut_x, cut_y (program-corrected values) to OUT structure
;2018-01-23	JD	- fixed a bug on #294, scoord undefined
;			- added correction of the last node location of the cut in piece-wise-linear cuts
;			- deleted everything that includes slit_x, slit_y, straight_slit, etc.
;2018-01-24	JD	- moved correction of the last cut node out of the if/endif for nodes GT 2
;			  (so that for cuts with 2 nodes the cut_x[1], cut_y[1] also gets corrected)
;
;2018-01-25	JD	- added /NAN keyword: if BOXCAR is set and there are NaN or Inf values, use /NAN
;			  note that running boxcar or nan will require significantly more time to complete
;			- added /EDGE_TRUNCATE to SMOOTH call
;2018-02-05	JD	- added handling of MIN_TIME and MAX_TIME with times only (does not require date)
;2018-02-06	JD	- corrected: min_time, max_time - anytim2utc conversion was missing
;			- added not(is_string(min_time & max_time)) to these keyword checks
;2018-02-22	JD	- added median(deriv(time.time))/2d (half-time resolution) to MIN_TIME, MAX_TIME checks
;			- Original MIN_TIME and MAX_TIME are now stored and returned (so that the input keyword is not overwritten)
;2018-07-23	JD	- removed the 'stop' command if both SLIT_WIDTH and SLIT_PIX are set
;			  (was producing an error in combination with stackplots.pro
;			   since the slit_pix was set in first stackplot then propagated to the oter 2)
;
;
; NOTES TO SELF
;
; after 2018-01-18	- slit_x, slit_y, straight_slit to be deprecated in the future (just get rid of them)


function ST_STACKPLOT, map, cut_x=cut_x, cut_y=cut_y, min_time=min_time, max_time=max_time, straight_cut=straight_cut,  $
; 			slit_x=slit_x, slit_y=slit_y, straight_slit=straight_slit,					$
			spline=spline, fast=fast, boxcar=boxcar, nan=nan, edge_truncate=edge_truncate, xp=xp, yp=yp,	$
			slit_width=slit_width, slit_pix=slit_pix, xap=xap, yap=yap

t0		= systime(/seconds)


;COORDINATE PROCESSING AND CHECKS
;-----------------------------------------------
get_map_coord, map[0], xm, ym
NX		= n_elements(map[0].data[*,0])
NY		= n_elements(map[0].data[0,*])
; NMAP		= n_elements(map[*].id)
DX		= map[0].dx
DY		= map[0].dy

; XM		= dblarr(NX,NY,NMAP)
; YM		= dblarr(NX,NY,NMAP)

; for im=0, NMAP-1, 1 do begin
;   get_map_coord, map[im], xmc, ymc
;   xm[*,*,im] = xmc
;   ym[*,*,im] = ymc
; endfor


; if keyword_set(slit_x) and not(keyword_set(cut_x)) then cut_x=slit_x
; if keyword_set(slit_y) and not(keyword_set(cut_y)) then cut_y=slit_y
; if keyword_set(straight_slit) and not(keyword_set(straight_cut)) then straight_cut=straight_slit

if keyword_set(XP) or keyword_set(YP) then begin
	print, ''
	print, '-!- ST_STACKPLOT: Keywords XP or YP not allowed to be used as inputs'
	print, ''
	stop
endif


if not(keyword_set(cut_X)) or not(keyword_set(cut_Y)) 	$
   or (n_elements(cut_X) NE n_elements(cut_Y))		then begin

	print, ''
	print, '-!- ST_STACKPLOT:  Malformed or missing cut_X or cut_Y coordinates'
	print, ''
	stop
endif

cut_X	 	= double(cut_X)
cut_Y		= double(cut_Y)

if (min(cut_X) LT min(xm)) or (max(cut_X) GT max(xm)) or $
   (min(cut_Y) LT min(ym)) or (max(cut_Y) GT max(ym)) then begin

	print, ''
	print, '-!- ST_STACKPLOT:  OUT of map ranges: Cut_X, Cut_Y'
	print, ''
	stop
endif


; if /straight_cut is set, treat the cut as if it has only 2 points
if keyword_set(straight_cut) then begin
  print, '% ST_STACKPLOT: Treating cut as if it has only 2 nodes '
  cut_X		= cut_X[0:1]
  cut_y		= cut_Y[0:1]
endif

if keyword_set(straight_cut) and (cut_X[0] GT cut_X[1]) then begin
	cut_X	= reverse(cut_X)
	cut_Y	= reverse(cut_Y)
	print, 'Reversing Cut_X and Cut_Y coordinates...'
endif
if keyword_set(straight_cut) and (cut_X[0] EQ cut_X[1]) then begin

	if (cut_Y[0] GT cut_Y[1]) then begin
		cut_Y	= reverse(cut_Y)
		print, 'Reversing Cut_Y coordinates...'
	endif
endif

; if (n_elements(cut_X) EQ 2) and (n_elements(cut_Y) EQ 2) then begin
;   straight_cut = 1
; endif


if keyword_set(slit_width) and keyword_set(slit_pix) then begin
	print, '-!- ST_STACKPLOT: Set only one of the SLIT_WIDTH or SLIT_PIX'
; 	stop
endif

if keyword_set(slit_pix) then begin
  if (fix(slit_pix) NE slit_pix) then begin
	slit_pix	= fix(slit_pix +0.5)
	print, '% ST_STACKPLOT: Keyword SLIT_PIX not an integer value. Corrected to '+trim(slit_pix)+' pixels'
  endif
endif



;TIME PROCESSING
;-----------------------------------------------

if not(keyword_set(min_time)) or not(is_string(min_time)) then begin
  print, '% ST_STACKPLOT:  MIN_TIME not set, assuming first map'
  min_time = map[0].time
endif
if not(keyword_set(max_time)) or not(is_string(min_time)) then begin
  print, '% ST_STACKPLOT:  MAX_TIME not set, assuming last map'
  max_time = map[n_elements(map.id)-1].time
endif

time		= anytim(map.time, out_style='vms')

; min_time_orig	= min_time
; max_time_orig	= max_time
min_time	= anytim(min_time, out_style='vms')
max_time	= anytim(max_time, out_style='vms')

if (strmid(min_time, 0, 11) EQ ' 1-Jan-1979') then strput, min_time, strmid(time[0],		      0, 11), 0
if (strmid(max_time, 0, 11) EQ ' 1-Jan-1979') then strput, max_time, strmid(time[n_elements(time)-1], 0, 11), 0

time		= anytim2utc(time)
min_time	= anytim2utc(min_time)
max_time	= anytim2utc(max_time)

;min_time & max_time OK?
if (min_time.time LT min(time.time)-median(deriv(time.time))/2d) or $
   (max_time.time GT max(time.time)+median(deriv(time.time))/2d) then begin
	print, ''
	if (min_time.time LT min(time.time)-median(deriv(time.time))/2d) then 	$
		print, '-!- ST_STACKPLOT:  INCORRECT time range: MIN_TIME'
	if (max_time.time GT max(time.time)+median(deriv(time.time))/2d) then	$
		print, '-!- ST_STACKPLOT:  INCORRECT time range: MAX_TIME'
	print, ''
	stop
endif

;midnight in the data?
; if (max_time.mjd EQ min_time.mjd +1L) then begin
; 	max_time.mjd	= min_time.mjd
; 	max_time.time	= max_time.time +3600L* 1000L* 24L
; endif
;
; ;more than one day?
; if (max_time.mjd GE min_time.mjd +2L) then begin
; 	print, ''
; 	print, '-!- ST_STACKPLOT:  INCORRECT time range: More than 2 dayap not currently implemented'
; 	print, ''
; 	stop
; endif

;midnight in the data?
Ndays	= (max_time.mjd - min_time.mjd)
if (Ndays GT 0) then begin
	max_time.mjd	= min_time.mjd
	max_time.time	= max_time.time +(3600L *1000L *24L *long(Ndays))
endif
if (Ndays GT 9) then begin
	print, ''
	print, '-!- ST_STACKPLOT: Too many days!', Ndays
	print, ''
	stop
endif

;correct for midnight in MAP times if neccessary
;untested!
if (Ndays GT 0) then begin
  for iday=1L, long(Ndays), 1L do begin
    inextday	= where(time.mjd EQ time[0].mjd +iday)
    time[inextday].mjd	= time[0].mjd
    time[inextday].time	= time[inextday].time +(3600L *1000L *24L *long(iday))
  endfor
endif

;get selected times
ind_t		= where((time.time/3.6d6 GE min_time.time/3.6d6) and $
			(time.time/3.6d6 LE max_time.time/3.6d6))
NTIMES		= n_elements(ind_t)



;STACKPLOTS - POINT IDENTIFICATION
;-----------------------------------------------

; has SPLINE been set?
if keyword_set(spline) then begin


	;curved cuts - with interpolation between points
	;do not know yet how to enforce strict 0.6 arcsec distance between points
	print, '% ST_STACKPLOT:  Creating spline curve through the cut nodes '
	spline_p, Cut_X, Cut_Y, xp, yp, interval=(DX+DY)/2d, /double
	npx		= n_elements(xp)
	DS		= 0d
	phi		= dblarr(npx)
	for i=1, npx-1 do begin
	  DS 		= DS +sqrt((xp[i]-xp[i-1])^2d +(yp[i]-yp[i-1])^2d)/double(npx-1d)
	  phi[i]	= atan( (yp[i]-yp[i-1])/(xp[i]-xp[i-1]) )
	endfor
	scoord		= dindgen(npx) *DS
	phi[0]		= phi[1]

endif else begin


	; Perform the point calculation for the first segment
	; this part of the code gets executed for both /STRAIGHT_CUT as well as for cuts with more than 2 nodes

	sgn	 	= +1d
	if (cut_X[0] NE cut_X[1]) then begin

	  ;non-vertical cuts
	  a		= (cut_Y[1] -cut_Y[0]) / (cut_X[1] -cut_X[0])
	  b		=  cut_Y[0] -a*cut_X[0]
	  dst		= sqrt( (cut_X[0] -cut_X[1])^2d +(cut_Y[0] -cut_Y[1])^2d )
	  DS		= 1d /sqrt( ((cos(atan(a)))/DX)^2d +((sin(atan(a)))/DY)^2d )
; 	  npx		= fix(dst /DS) +1
	  npx		= fix(dst /DS +0.5d) +1
	  phi		= replicate(atan(a), npx)

          if (cut_Y[0] EQ cut_Y[1]) and (n_elements(cut_x) EQ 2) then begin
	    scoord 	= dindgen(npx) *DS +cut_X[0]
	  endif else begin
	    scoord	= dindgen(npx) *DS
	  endelse

	  if (cut_x[1] LT cut_x[0]) then sgn=-1d
	  xp		= sgn *DX *dindgen(npx) *cos(phi[0]) +cut_X[0]			; if tan(phi[0]) = a, phi = <0,!dpi/2)
	  yp		= sgn *DY *dindgen(npx) *sin(phi[0]) +cut_Y[0]
	endif else begin

	  ;vertical cuts
; 	  npx		= fix(abs(Cut_Y[1] -Cut_Y[0]) /DY +1e-3) +1
	  npx		= fix(abs(Cut_Y[1] -Cut_Y[0]) /DY +0.5d) +1
	  if (cut_y[1] LT cut_Y[0]) then sgn = -1d
	  xp		= replicate(cut_X[0], npx)
	  yp		= sgn *DY *dindgen(npx) +cut_Y[0]

	  if keyword_set(straight_cut) or (n_elements(cut_y) EQ 2) then begin
	    scoord	= DY *dindgen(npx) +cut_Y[0]
	  endif else begin
	    scoord	= DY *dindgen(npx)
	  endelse

	  DS		= DY
  	  phi		= replicate(sgn *!dpi/2d, npx)
	endelse

; 	;correct cut_x[1], cut_y[1]
; 	if (cut_X[1] NE xp[n_elements(xp)-1]) or (cut_Y[1] NE yp[n_elements(yp)-1]) then begin
; 	  print, '% ST_STACKPLOT: Correcting cut_X[1]   from ' $
; 		+trim(cut_X[1])+' to '+trim(xp[n_elements(xp)-1])
; 	  print, '% ST_STACKPLOT: Correcting cut_Y[1]   from ' $
; 		+trim(cut_y[1])+' to '+trim(yp[n_elements(yp)-1])
; 	  cut_X[1]	= xp[n_elements(xp)-1]
; 	  cut_Y[1]	= yp[n_elements(yp)-1]
;          endif


	; does the cut have more than 2 points?
	; if yes, perform the analogous segment calculation for the next segments
	; note this part of the code does not get executed if /straight_cut is set
	if (n_elements(cut_x) GT 2) and (n_elements(cut_y) GT 2) then begin

	  for isegment=2, n_elements(cut_x)-1, 1 do begin

	    print, '% ST_STACKPLOT: Treating segment ' +trim(isegment)

	    ;use the last xp and yp as the start for the next segment
	    ;(this modifies the cut node location by less than 1 pixel)
	    if (cut_X[isegment-1] NE xp[n_elements(xp)-1]) or (cut_Y[isegment-1] NE yp[n_elements(yp)-1]) then begin
	      print, '% ST_STACKPLOT: Correcting cut_X['+trim(isegment-1)+']   from ' $
			+trim(cut_X[isegment-1])+' to '+trim(xp[n_elements(xp)-1])
	      print, '% ST_STACKPLOT: Correcting cut_Y['+trim(isegment-1)+']   from ' $
			+trim(cut_Y[isegment-1])+' to '+trim(yp[n_elements(yp)-1])

	      cut_X[isegment-1]	= xp[n_elements(xp)-1]
	      cut_Y[isegment-1]	= yp[n_elements(yp)-1]
            endif

	    sgn 	= +1d

	    ;now proceed analogously to the first segment
	    ;but keep the DS from the first segment for consistency
	    if (cut_X[isegment-1] NE cut_X[isegment]) then begin

		if (cut_x[isegment] LT cut_x[isegment-1]) then sgn=-1d

		a	= (cut_Y[isegment] -cut_Y[isegment-1]) / (cut_X[isegment] -cut_X[isegment-1])
		b	=  cut_Y[isegment-1] -a*cut_X[isegment-1]
		dst	= sqrt( (cut_X[isegment-1] -cut_X[isegment])^2d +(cut_Y[isegment-1] -cut_Y[isegment])^2d )
; 	        DS	=  1d /sqrt( ((cos(atan(a))/DX)^2d +((sin(atan(a))/DY)^2d )
		npx2	= fix(dst /DS +0.5d) +1
		scoord	= [scoord, 		max(scoord) +(dindgen(npx2-1)+1d) *DS]
		phi	= [phi,			replicate(atan(a), npx2-1)]

		if (cut_Y[isegment-1] EQ cut_Y[isegment]) then begin
		  scoord= scoord +cut_X[isegment-1]
		endif
		xp	= [xp, 			sgn *DX *(dindgen(npx2-1)+1d) *cos(phi[npx]) +cut_X[isegment-1] ]
		yp	= [yp, 			sgn *DY *(dindgen(npx2-1)+1d) *sin(phi[npx]) +cut_Y[isegment-1] ]
	    endif else begin

	    ;vertical cuts
		if (cut_Y[isegment] LT cut_y[isegment-1]) then sgn =-1d

		npx2	= fix(abs(Cut_Y[isegment] -Cut_Y[isegment-1]) /DS +0.5d) +1
		xp	= [xp, 			replicate(cut_X[isegment-1], npx2-1)]
		yp	= [yp, 			sgn *DS *(dindgen(npx2-1)+1d) +cut_Y[isegment-1] ]
		scoord	= [scoord, 		max(scoord) +DS *(dindgen(npx2-1)+1d) ]
		phi	= [phi, 		replicate(sgn *!dpi/2d, npx2-1) ]
	    endelse
	    npx		= npx +npx2 -1
  	  endfor

	endif else begin
	  isegment = 2
	endelse

	; finally, correct the last node of the cut
	if (cut_X[isegment-1] NE xp[n_elements(xp)-1]) or (cut_Y[isegment-1] NE yp[n_elements(yp)-1]) then begin
	  print, '% ST_STACKPLOT: Correcting cut_X['+trim(isegment-1)+']   from ' $
		+trim(cut_X[isegment-1])+' to '+trim(xp[n_elements(xp)-1])
	  print, '% ST_STACKPLOT: Correcting cut_Y['+trim(isegment-1)+']   from ' $
		+trim(cut_Y[isegment-1])+' to '+trim(yp[n_elements(yp)-1])
	  cut_X[isegment-1]	= xp[n_elements(xp)-1]
	  cut_Y[isegment-1]	= yp[n_elements(yp)-1]
	endif
endelse


; Is there finite a slit width?
; if yes, then define the slit (aperture)
; and calculate everything in pixels across the cut
if keyword_set(slit_width) then slit_pix = fix(slit_width/DS +0.5)

if keyword_set(slit_pix) then begin

  slit_width	= slit_pix *DS
  print, '% ST_STACKPLOT: Slit_width corrected to '+trim(slit_width)+' arc sec'

  xap		= dblarr(npx,2*slit_pix+1)
  yap		= dblarr(npx,2*slit_pix+1)

  for i=0, npx-1, 1 do begin
   for j=-slit_pix, +slit_pix, 1 do begin
	xap[i,j+slit_pix]	= xp[i] +j*DS*sin(phi[i])
	yap[i,j+slit_pix]	= yp[i] -j*DS*cos(phi[i])					; minus?
    endfor
  endfor
endif else begin
  xap		= xp								; ?
  yap		= yp
endelse

; for output
if NOT(keyword_set(slit_pix)) then slit_width = 0d

NSLIT_PIX	= n_elements(xap[0,*])
stack_slit	= dblarr(npx,NSLIT_PIX,NTIMES)



;STACKPLOTS
;-----------------------------------------------

dt			= median(deriv(time.time))/1d3
out			= { stackplot:dblarr(npx,NTIMES),   id:map[0].id, tcoord:map[ind_t].time, dt:dt, 		$
			    scoord:scoord, ds:ds, xp:xp, yp:yp, cut_x:cut_x, cut_y:cut_y, spline:keyword_set(spline), 	$
			    xap:xap, yap:yap, slit_width:slit_width}

for it=0, NTIMES-1, 1 do begin
  if NOT(keyword_set(FAST)) then get_map_coord, map[ind_t[it]], xm, ym
  if not(keyword_set(boxcar)) then begin
;     out.stackplot[*,it] = interpolate(map[ind_t[it]].data, $
; 				      interpol(dindgen(NX), xm[*,0],xp), $
; 				      interpol(dindgen(NY), ym[0,*],yp) )
    for j=0, NSLIT_PIX-1 do begin
      stack_slit[*,j,it]= interpolate(map[ind_t[it]].data, $
				      interpol(dindgen(NX), xm[*,0], xap[*,j]), $
				      interpol(dindgen(NY), ym[0,*], yap[*,j]) )
    endfor
  endif else begin
;     out.stackplot[*,it] = interpolate(smooth(map[ind_t[it]].data, boxcar), $
; 				      interpol(dindgen(NX), xm[*,0],xp), $
; 				      interpol(dindgen(NY), ym[0,*],yp) )
    for j=0, NSLIT_PIX-1 do begin
      stack_slit[*,j,it]= interpolate(smooth(map[ind_t[it]].data, boxcar, nan=nan, edge_truncate=edge_truncate), $
				      interpol(dindgen(NX), xm[*,0], xap[*,j]), $
				      interpol(dindgen(NY), ym[0,*], yap[*,j]) )
    endfor
  endelse
;   out.stackplot[*,it] = total(stack_slit[*,*,it], 2, /double)/double(NSLIT_PIX)
endfor

out.stackplot 	= total(stack_slit, 2, /double)/double(NSLIT_PIX)


; min_time	= min_time_orig
; max_time	= max_time_orig
min_time	= anytim(min_time, out_style='vms')
max_time	= anytim(max_time, out_style='vms')


print, '% ST_STACKPLOT: Took ' +string((systime(/seconds)-t0), format='(F7.3)') +' seconds to complete.'
return, out
end
