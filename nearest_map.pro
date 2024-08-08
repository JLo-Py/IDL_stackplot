function nearest_map, maps_in, ref_time=ref_time, ref_map=ref_map, quiet=quiet, earlier=earlier, later=later

;INPUT
;maps_in			- maps for which nearest neighbor to the reference map is to be searched

;ref_time			- reference time.
;				  Has to be in a format compatible with map time tag; e.g. '2014-10-09 16:40:00'
;OR
;ref_map			- reference map
;
;OUTPUT
;ind				- index of the maps_in closest to the time of ref_map
;
;
;HISTORY
;2014-12-03	JD		- adapted from nearest_time_neighbor
;2015-02-26	JD		- added keyword QUIET
;				- added keywords EARLIER and LATER: only to find the nearest time from earlier/later than the ref_time
;
;2017-10-11	JD		- renamed to NEAREST_MAP.PRO
;2018-02-21	JD		- changed print output to include the function name
;2018-08-30	JD		- added handling of REF_TIME as in st_stackplot.pro: 
;				  with times only (does not require date)
;2019-02-20	JD		- backed up
;				  CHANGED so as to handle arrays of ref_time of ref_maps
; 
;2020-03-13	JD		- [coronavirus homeoffice]
;				  JL discovered a bug: The program cannot discern between different dates present in the input
;				  (not surprising, since I never really tried)
;				  --> copy-pasted & adapted a piece of code from st_stackplot.pro
;				- also adapted the date_match (now obsolete?) for possible use with IDL 8.x



if NOT(keyword_set(ref_time)) and NOT(keyword_set(ref_map)) then begin
  print, '-!- NEAREST_MAP:  No reference time or map set'
  stop
endif
if keyword_set(ref_time) and keyword_set(ref_map) then begin
  print, '-!- NEAREST_MAP:  Both reference time and reference map set'
  stop
endif

if keyword_set(ref_map)  then ref_time = anytim2utc(ref_map.time)
if keyword_set(ref_time) then ref_time = anytim2utc(ref_time)


FOR it=0L, n_elements(ref_time)-1L, 1L DO BEGIN

  reft		= anytim(ref_time[it], out_style='vms')
  if (strmid(reft, 0, 11) EQ ' 1-Jan-1979') then strput, reft, strmid(maps_in[0].time,  0, 11), 0
  if keyword_set(ref_time) then reft = anytim2utc(reft)


  N_MAPS		= n_elements(maps_in[*].id)
  ttcoord		= anytim2utc(maps_in[*].time)
  min_time		= anytim2utc(min(maps_in[*].time))
  max_time		= anytim2utc(max(maps_in[*].time))


  if (float(!version.release) GE 8.) then begin
    date_match		= where(ttcoord[*].mjd eq reft.mjd, /null)
    if (date_match EQ !null) then print, '% NEAREST_MAP:  Maps from different dates'
  endif else begin
    date_match		= where(ttcoord[*].mjd eq reft.mjd)
    if (date_match[0] eq -1) then  print, '% NEAREST_MAP:  Maps from different dates'
  endelse


  ; correct for different dates if necessary
  ; adapted from st_stackplot and plot_st_stackplot (still untested!)
  Ndays		= (max_time.mjd - min_time.mjd)
  if (Ndays GT 0) then begin
	max_time.mjd	= min_time.mjd
	max_time.time	= max_time.time +(3600L *1000L *24L *long(Ndays))
  endif
  if (Ndays GT 9) then begin
	print, ''
	print, '-!-  NEAREST_MAP:  Too many days!', Ndays
	print, ''
	stop
  endif
  if (Ndays GT 0) then begin
    for iday=1L, long(Ndays), 1L do begin
	inextday		= where(ttcoord.mjd EQ ttcoord[0].mjd +iday)
	ttcoord[inextday].mjd	= ttcoord[0].mjd
	ttcoord[inextday].time	= ttcoord[inextday].time +(3600L *1000L *24L *long(iday))
    endfor
  endif

  ; now correct the reference time
  if (reft.mjd NE min_time.mjd) then begin
    Ndays_reft	= reft.mjd -min_time.mjd
    print, '% NEAREST_MAP: The reference date is  '+trim(Ndays_reft)+' days later than the first map'
    reft.time	= reft.time +(3600L *1000L *24L *long(Ndays_reft))
  endif



  if keyword_set(earlier) 	then ind = where(ttcoord[*].time LE reft.time)
  if keyword_set(later)   	then ind = where(ttcoord[*].time GE reft.time)
  if not(keyword_set(earlier)) and $
     not(keyword_set(later))	then ind = indgen(N_MAPS)


  dumb		= min(abs(ttcoord[ind].time - reft[0].time), C)
  sgn		= dumb/(ttcoord[ind[C]].time - reft[0].time)

  if (it EQ 0) then indt = ind[C] else indt = [indt, ind[C]]

  if NOT(keyword_set(quiet)) then print, '% NEAREST_MAP:  Minimum time difference found is '+trim(sgn*dumb/1d3)+' seconds'

ENDFOR

return, indt
END
