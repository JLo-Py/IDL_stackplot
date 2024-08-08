PRO plot_cut_spline, st_stackplot, thick=thick, hsize=hsize, charsize=charsize, charthick=charthick, color=color, len_tick=len_tick, interval=interval, no_arrow=no_arrow, srange=srange, notickmarks=notickmarks, text_o=text_o, dev_fact=dev_fact, shx_txt=shx_txt, shy_txt=shy_txt, change_orien=change_orien, scraplast=scraplast

;PURPOSE
;plots a cut onto a map
;
;INPUT
;st_stackplot		- stackplot created by st_stackplot with keyword spline = 1
;
;OPTIONAL INPUT
;Color			- color # used to plot the cut. Default = 255
;Thick			- arrow thick size
;Hsize			- arrow head size
;Charsize, Charthick	- character size & thickness
;Text_o                 - orientation of the captions 

;tsx and tsy            - text shifts, still need to be coded!

;HISTORY
;2020-12-24	JL	- written, based on JD's plot_cut.pro routine
;2020-12-25     JL      - program now handles plotting of major tickmarks
;2020-12-26     JL      - in both X and PS devices, solid lines 'oplotted' over images are ugly - perhaps a plotting or rendering issue???
;                       - ticks are now plotted via the plots function - BORING
;2021-01-02     JL      - workaround for a bug concerning the caption of the cut at the edge node

if not(keyword_set(st_stackplot.spline)) then begin
  message, '-!- PlOT_CUT_SPLINE now adapted for splined cuts only'
endif

if keyword_set(srange) and n_elements(srange) NE 2 then begin

  print, ' -!- keyword SRANGE not properly formatted:', srange
  print, '     Should be 2-element vector'
  stop
  
endif

if not(keyword_set(srange)) then begin
  srange = [min(st_stackplot.scoord), max(st_stackplot.scoord)]
endif

if not(keyword_set(charthick))	then charthick=2
if not(keyword_set(color))	then color=0
if not(keyword_set(thick))	then thick=3

if (!d.name EQ 'X') or (!d.name EQ 'WIN') then begin

    if not(keyword_set(hsize))	    then hsize=30
    if not(keyword_set(charsize))   then charsize=2.5
    if not(keyword_set(dev_fact))   then dev_fact=1.
    dev_tick    = 150

endif else begin

    if not(keyword_set(hsize))	    then hsize=300
    if not(keyword_set(charsize))   then charsize=1.25
    if not(keyword_set(dev_fact))   then dev_fact=10.
    dev_tick    = 100

endelse

if not(keyword_set(len_tick)) 	then len_tick	= sqrt((!x.crange[1]-!x.crange[0])^2 +(!y.crange[1]-!y.crange[0])^2d) /dev_tick

;1. extract the splined cut

inds		= where((st_stackplot.scoord GE srange[0]) and (st_stackplot.scoord LE srange[1]))
scoord		= st_stackplot.scoord[inds]

xp		= st_stackplot.xp[inds]
yp		= st_stackplot.yp[inds]

;2. set tick interval (if not set manually)

ncoord		= n_elements(scoord)
cut_length	= scoord[ncoord-1]-scoord[0]		; in arc sec
precision	= 1d-2					; precision of the point; [arc sec]

if not(keyword_set(interval)) then begin 

  if (cut_length LE 1000d) then interval = 100 else interval = 200
  if (cut_length LE 250d)  then interval = 50
  if (cut_length LE 25d)   then interval = 10
  if (cut_length LE 10d)   then interval = 3
  
endif

;3. interpolate the cut to a finer grid, where deriv(cut_interpol) < precision

steps=deriv(scoord)
old_grid=steps[0]/precision
new_grid=fix(ncoord*old_grid)

cl_int=interpol(scoord, new_grid)
xp_int=interpol(xp, new_grid)
yp_int=interpol(yp, new_grid)

;4. get the interval points 

if keyword_set(scraplast) then begin 
    np=round(cut_length/interval)
    endif else begin 
    np=round(cut_length/interval)+1          ;including the tick underneath the arrow 
endelse 

if srange[0] NE 0. then begin 
    np=round(cut_length/interval)
    message, '!!!   TO BE CODED! Node calculation does not work with srange set   !!!'
endif 

clr=round(cut_length)+1                 ;round of length
cut_int=fix(interpol([0,clr], clr+1))   ;integers of length

nodes=intarr(np)                        ;integers, make them strings later on

; if not(keyword_set(tick_o)) then begin 
;     print, '!!!   Tick orientation automatically set to 1,  EXPECT STRANGE RESULTS   !!!'
;     tick_o=intarr[np]
;     tick_o[*]=1
; endif  

if (np-1)*interval GT clr then begin 
    for ip=0, np-2 do nodes[ip]=cut_int[ip*interval]        ;filling the last element -- might not be needed once /scraplast is set?
endif else begin 
    for ip=0, np-1 do nodes[ip]=cut_int[ip*interval]
endelse 
    
;5. find the integers corresponding to your interval's points in the interpolated array

point_ind=intarr(np)

for ipi=0, np-1 do point_ind[ipi]=where((cl_int GE nodes[ipi]-precision/2) and (cl_int LE nodes[ipi]+precision/2))

;6. and their X and Y positions 

xpos_n= xp_int[point_ind] 
ypos_n= yp_int[point_ind] 

;7. MAIN PLOTTING: plot the line ended with an arrow... 
;does not need to be in a cycle for now... but watch for the last interval, so it does not look like in JD's plotcut 

plots, xp, yp, thick=thick, color=color, /data, linestyle=0 
arrow, xp[ncoord-2], yp[ncoord-2], xp[ncoord-1], yp[ncoord-1], thick=thick, color=color, /data, /solid, hsize=hsize

;8. plot the MAJOR ticks and coordinate for each segment individually

nodes_str   =string(nodes)
nodes_str   =trim(nodes_str)

for ti=0, np-1 do begin 

    if point_ind[ti] EQ 0 then begin 
            ind_lo  = 0
            ind_up  = 2
            point_ind[ti] = 1                   ;if node is at scoord=0, shift it by one pixel in the interpolated grid
    endif else begin 
            ind_lo=point_ind[ti]-1              ;points closest to the point - your tangent will cross these 
            ind_up=point_ind[ti]+1
    endelse 
    
    line=linfit([xp_int[ind_lo],xp_int[ind_up]], [yp_int[ind_lo],yp_int[ind_up]])
    
    minx=min(st_stackplot.xp)
    maxx=max(st_stackplot.xp)*2                 ;widen the range of x in case one needs longer tickmarks! 
    disx=abs(maxx)-abs(minx)
    
    tickgrid=fix(disx/precision)
    
    tickgrid=float(tickgrid)*2                  ;the grid of x and then y_p must be very fine, so the tick 'touches' the cut

    x=interpol([minx, maxx], tickgrid)      
    
    y_l=line[1]*x+line[0]                       ;equation of a line tangent to the node
    
    point_p=[xpos_n[ti],ypos_n[ti]]             ;produce a line perpendicular to this one crossing the node

    y_p=line_perp(x, line[1], line[0], /perp, point=point_p, m=m, n=n)

    orientation	= atan(m)*180/!dpi
    
    print, 'Tick orientation is: '+trim(string(orientation))+''
    sign = (-1d)*(orientation LT 0d) + (+1d)*(orientation GE 0d)
    
    if sign EQ 1 then begin 
        id=where(y_p GE point_p[1])         ;0.3 so ticks will reach into the cut!
    endif else begin 
        id=where(y_p LE point_p[1])
    endelse 
        
    xid=x[id]
    yid=y_p[id]

    dists=fltarr(n_elements(id))

    for indd=0L, n_elements(id)-1 do begin      
        data=[[point_p],[xid[indd],yid[indd]]]     
        dists[indd]=distance_measure(data)     
    endfor 
    
    differences=abs(dists-len_tick*dev_fact)
    d=where(min(differences) EQ differences)
    
    if keyword_set(change_orien) then begin 
    
        message, 'To be coded! '
    
    endif else begin 
    
;         oplot, xid[0:d], yid[0:d], thick=thick, color=color
        plots, [point_p[0], xid[d]], [point_p[1], yid[d]], thick=thick, color=color, /data, linestyle=0 
    
    endelse 

    if keyword_set(shx_txt) then begin 
    
        if (n_elements(shx_txt) NE np) or (n_elements(shy_txt) NE np) then message, 'ERROR:    Number of nodes does not equal the number of text positions.'
    
        xyouts, xid[d]+shx_txt[d], yid[d]+shy_txt[d], /data, nodes_str[ti], color=color, charsize=charsize, charthick=charthick, orientation=text_o
    
    endif else begin 
        
        shiftx_txt=sign*len_tick*dev_fact*sin(orientation)
        shifty_txt=sign*len_tick*dev_fact*cos(orientation)
        
        xyouts, xid[d]+shiftx_txt, yid[d]+shifty_txt, /data, nodes_str[ti], color=color, charsize=charsize, charthick=charthick, orientation=text_o

    endelse 
endfor 

end 
