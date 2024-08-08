PRO PLOT_CUT, st_stackplot, thick=thick, hsize=hsize, charsize=charsize, charthick=charthick, color=color, text_color=text_color, text_out=text_out, text_x=text_x, text_y=text_y, orientation=orientation, len_tick=len_tick, interval=interval, no_arrow=no_arrow, srange=srange, notickmarks=notickmarks, dev_fact=dev_fact, nominor=nominor


;PURPOSE
;plots a cut onto a map
;
;INPUT
;st_stackplot		- stackplot created by st_stackplot
;
;OPTIONAL INPUT
;TEXT			- text to write
;Text_x, Text_y		- lower left coordinated for the text
;Color			- color # used to plot the cut. Default = 255
;Thick			- arrow thick size
;Hsize			- arrow head size
;Charsize, Charthick	- character size & thickness
;
;NOTE
; This program is ONLY to be used on the same window / eps file as with plotmap
; (otherwise expect strange results)
;
;OPTIONAL OUTPUT
;Orientation		- orientation for the text
;
;
;HISTORY
;2015-12-16	JD	- written (adapted from a piece of code in flares/pro/plot_2014-09-10_aia.pro for 171A
;2015-01-11	JD	- added keyword NO_ARROW. Default is to plot the cut as an arrow
;			- renamed keyword TEXT to TEXT_OUT
;			- changed the default color to 0 (otherwise setting color=0 would NOT work)
; 			- changed len_tick modifier to 100 (was 150)
;			- added keyword TEXT_COLOR
;2016-01-14	JD	- changed len_tick to !p.ticklen instead of calculating it from x/y.crange
;			  reverted back, does not work for PS (?)
;			- changed how dev_tick is computed: THIS WILL REQUIRE SOME REDESIGN IN THE FUTURE
;			- added keyword SRANGE
;			- added keyword NOTICKMARKS: only cut & ticks will be plotted
;2016-01-18	JD	- modified the condition for srange: has to be keyword_set(srange) and n_elements(srange) NE 2
;			  added automatic setting of srange if not set by keyword
;2016-01-19	JD	- added DEV_FACT as a keyword
;2016-11-16	JD	- added INTERVAL as a keyword
;2017-05-10	JD	- changed the variable 'ia' to longword type
;2018-01-23	JD	- moved the interval calculation
;			- attempted rewrite to plot the piecewise-linear cuts
;2018-01-24	JD	- finished
;			- added checks for splined cuts (this program is not yet ready to plot these cuts)
;			- corrected a bug with NNODE and cut_x, cut_y calculation: now uses cut_xt, cut_yt

;2020-12-25     JL      - added keyword nominor - if set, minor tickmarks are not plotted 
       
if keyword_set(st_stackplot.spline) then begin
  print, '-!- PLOT_CUT: Plot_cut not yet adapted to /SPLINE cuts'
  stop
endif

if keyword_set(srange) and n_elements(srange) NE 2 then begin
  print, ' -!- keyword SRANGE not properly formatted:', srange
  print, '     Should be 2-element vector'
  stop
endif
if not(keyword_set(srange)) then begin
  srange = [min(st_stackplot.scoord), max(st_stackplot.scoord)]
  ; note that this is different than in the plot_st_stackplot, where srange refers to the srange for the plotting of the stackplot
  ; this corresponds to the center of the pixels only
endif


if not(keyword_set(charsize))	then charsize=!p.charsize
if not(keyword_set(charthick))	then charthick=!p.charthick

if (!d.name EQ 'X') or (!d.name EQ 'WIN') then begin

  if not(keyword_set(color))	then color=0
  if not(keyword_set(thick))	then thick=!p.thick
  if not(keyword_set(hsize))	then hsize=10

  dev_tick		= 150
  if not(keyword_set(dev_fact)) then dev_fact= 2.5
endif else begin

  if not(keyword_set(color))	then color=0
  if not(keyword_set(thick))	then thick=4
  if not(keyword_set(hsize))	then hsize=300

  dev_tick		= 100
  ;   dev_fact		= 3.0
  if not(keyword_set(dev_fact)) then dev_fact= 2.0 *sqrt((!x.crange[1]-!x.crange[0])^2d +(!y.crange[1]-!y.crange[0])^2d)/100d < 5.0
;   print, dev_fact
endelse


if not(keyword_set(text_color)) then text_color = color
if not(keyword_set(len_tick)) 	then len_tick	= sqrt((!x.crange[1]-!x.crange[0])^2d 		$
						      +(!y.crange[1]-!y.crange[0])^2d) /dev_tick
; if not(keyword_set(len_tick)) 	then len_tick	= !p.ticklen
; print, len_tick
;this would not work if the plotmap is not set (if no plotting window exists)


; ;set device factor for the distance between the tick and the tick name
; if (!d.name EQ 'X') or (!d.name EQ 'WIN') then 	dev_fact = 1.5
; if (!d.name EQ 'PS') then 			dev_fact = 3.0






; extract the cut
inds		= where((st_stackplot.scoord GE srange[0]) and (st_stackplot.scoord LE srange[1]))
scoord		= st_stackplot.scoord[inds]
xp		= st_stackplot.xp[inds]
yp		= st_stackplot.yp[inds]

cut_x		= st_stackplot.cut_x
cut_y		= st_stackplot.cut_y

NCUT		= n_elements(cut_x)
NNODE		= 0


; set tick interval (if not set already)
NSCOORD		= n_elements(scoord)
cut_length	= scoord[NSCOORD-1]-scoord[0]		; in arc sec
precision	= 1d-2					; precision of the point; [arc sec]
; NP		= long(cut_length/precision)

if not(keyword_set(interval)) then begin 
  if (cut_length LE 1000d) then interval = 100 else interval = 200
  if (cut_length LE 250d)  then interval = 50
  ; if (cut_length LE 100d)  then interval = 10
  if (cut_length LE 25d)   then interval = 10
  if (cut_length LE 10d)   then interval = 1		; unlikely to be this short, but still keep it
endif

;larger are not necessary - R_Sun
sub_interval	= interval / 5





; fit what cut points are inside [srange[0], srange[1]]
for icut=0L, NCUT-1L, 1L do begin

  if (float(!version.release) GE 8.) then begin
    dummy	= where( ((xp EQ cut_x[icut]) and (yp EQ cut_y[icut])), count, /null)
  endif else begin
    dummy	= where( ((xp EQ cut_x[icut]) and (yp EQ cut_y[icut])), count)
  endelse

;   if (count GT 0) and (count NE 1) then begin
;     print, '-!- PLOT_CUT: Something is rotten in the state of Denmark. Make sure to tidy up.'
;     stop
;   endif

  NNODE		= NNODE +count
  
  ; if the first node is reached, cull the (cut_x, cut_y) from the left 
  if (NNODE EQ 1) and (count NE 0) then begin
    cut_xt	= cut_x[icut:NCUT-1]
    cut_yt	= cut_y[icut:NCUT-1]
  endif
endfor

; if there is less than 1 node point inside, keep the original
; and for all intents and purposes treat it as a straight cut with 2 nodes
if (NNODE GE 1) then begin
  cut_xt	= cut_xt[0:NNODE-1]
  cut_yt	= cut_yt[0:NNODE-1]
  if (cut_xt[0] NE xp[0]) and (cut_yt[0] NE yp[0]) then begin
    cut_xt	= [xp[0], cut_xt]
    cut_yt	= [yp[0], cut_yt]
  endif
  if (cut_xt[n_elements(cut_xt)-1] NE xp[n_elements(xp)-1]) and $
     (cut_yt[n_elements(cut_yt)-1] NE yp[n_elements(yp)-1]) then begin
    cut_xt	= [cut_xt, xp[n_elements(xp)-1]]
    cut_yt	= [cut_yt, yp[n_elements(yp)-1]]
  endif
endif else begin
  cut_xt	= [xp[0], xp[n_elements(xp)-1]]
  cut_yt	= [yp[0], yp[n_elements(yp)-1]]
endelse
NNODE		= n_elements(cut_xt)

; print, NNODE, cut_xt, cut_yt

; now plot the cut piece-wise
for inode=1, NNODE-1 do begin

  if (float(!version.release) GE 8.) then begin
    inode_a	= where( ((xp EQ cut_xt[inode-1]) and (yp EQ cut_yt[inode-1])), /null )
    inode_b	= where( ((xp EQ cut_xt[inode])   and (yp EQ cut_yt[inode]  )), /null )
  endif else begin
    inode_a	= where( ((xp EQ cut_xt[inode-1]) and (yp EQ cut_yt[inode-1])) )
    inode_b	= where( ((xp EQ cut_xt[inode])   and (yp EQ cut_yt[inode]  )) )
  endelse

  NS		= n_elements(scoord[inode_a:inode_b])

  a 		= linfit(xp[inode_a:inode_b], yp[inode_a:inode_b])
  orientation	= atan(a[1])*180/!dpi
  sign		= (-1d)*(orientation LT 0d) + (+1d)*(orientation GE 0d)

  piece_length	= (scoord[inode_b] -scoord[inode_a])[0]				; in arc sec
  NP		= long(piece_length/precision)

 
  ;plot the cut
  
 print, XP[inode_b], YP[inode_b]
 
 u=convert_coord([!D.X_SIZE,!D.X_SIZE-hsize], /device, /to_data)    ;you have to treat margins
 
 print, u
;  
;  arx=XP[inode_b]
;  ary=YP[inode_b]
;  
;  sx=
;  
  if (keyword_set(no_arrow)) or (inode NE NNODE-1) then begin
    plots, [XP[inode_a], XP[inode_b]], [YP[inode_a], YP[inode_b]], /data, thick=thick, color=color
  endif else begin
    arrow,  XP[inode_a], YP[inode_a], XP[inode_b], YP[inode_b], /data, /solid, thick=thick, color=color, hsize=hsize
; ;     arrow,  XP[inode_a], YP[inode_a], XP[inode_b+hsize/2.], YP[inode_b+hsize/2.], /data, /solid, thick=thick, color=color, hsize=hsize


  endelse

  ;determine cut orientation
  if (cut_x[inode] LT cut_x[inode-1]) then begin
    sgn 	= -1d
    xtx		= xp[inode_b]
    ytx		= yp[inode_b]
  endif else begin
    sgn 	= +1d
    xtx		= xp[inode_a]
    ytx		= yp[inode_a]
  endelse
    
  ;plot the text if required
  if keyword_set(text_out) then begin
  ;   len_win	= sqrt(!x.crange[1]-!x.crange[0])^2d +(!y.crange[1]-!y.crange[0])^2d)
    if not(keyword_set(text_x)) and (inode EQ 1) then text_x = xtx +!d.x_size*1d-3
    if not(keyword_set(text_y)) and (inode EQ 1) then text_y = ytx +!d.x_size*1d-3
    if (inode EQ 1) then xyouts, text_x, text_y, /data, text_out, orientation=orientation, $
				 color=text_color, charsize=charsize, charthick=charthick
  endif
  
  ;plot axis along the cut
  nprec		= strlen(trim(precision))
  scoord_a	= double( string(scoord[inode_a], format='(F10.'+trim(nprec-2)+')') )
;   scoord_a	= double( string(scoord[inode_a], format='(F'+trim(nprec)+'.'+trim(nprec-2)+')') )
;   scoord_a	= double(trim( scoord[inode_a] ))
    if ((scoord_a mod sub_interval) EQ 0d) then step = sub_interval else step = 1
    
    if keyword_set(nominor) then minfact=1e3 else minfact=2.5
  
;   for ia=0L, NP-1L, sub_interval do begin
  for ia=0d, NP-1d, step do begin
    ;ticks
    if (((ia*precision +scoord_a) mod interval) EQ 0d) then length=len_tick*charsize else length = len_tick*charsize/minfact
    if (((ia*precision +scoord_a) mod sub_interval) EQ 0d) then begin
       	 plots, (XP[inode_a] +sgn *ia*precision*cos(atan(a[1])))[0] +[0, +length*sin(atan(a[1]))],			$
		(YP[inode_a] +sgn *ia*precision*sin(atan(a[1])))[0] +[0, -length*cos(atan(a[1]))],			$
	   	/data, thick=thick, color=color
    endif
  
;     print, inode, piece_length, NP, ia, sub_interval, interval, ia*precision, scoord[inode_a], scoord_a, ia*precision +scoord_a, $
; 	   ((ia*precision +scoord_a) mod interval), ((ia*precision +scoord_a) mod sub_interval)
  
    ;tickmarks
    if not(keyword_set(notickmarks)) then begin 
      if (((ia*precision +scoord_a) mod interval) EQ 0d) then begin
  	    xyouts, XP[inode_a] +sgn*(ia*precision*cos(atan(a[1]))						$
					  +sgn*length*(dev_fact*sin(atan(a[1])) -cos(atan(a[1])))),		$
		    YP[inode_a] +sgn*(ia*precision*sin(atan(a[1]))						$
					  -sgn*length*(dev_fact*cos(atan(a[1])) +sin(atan(a[1])))),  /data,	$
		    trim(ia*precision +scoord_a), orientation=orientation,					$
		    charsize=charsize, charthick=charthick, color=color
      endif
    endif
  endfor

endfor



END
