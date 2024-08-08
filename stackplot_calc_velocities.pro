PRO stackplot_calc_velocities, line_x=line_x, line_y=line_y, quiet=quiet, transpose=transpose, ds=ds, dt=dt, soho=soho, start_time=start_time, vel=vel, d_vel=d_vel

;2016-01-11	JD	- written as a quick utility to get the velocities & uncertainties from measurements done previously
;2016-01-19	JD	- updated to reflect the latest stackplot_cursor.pro
;2016-01-20	JD	- changed the calculation of D_vel using the sqrt sigma formula
;2016-02-22	JD	- added keyword start_time
;			- changed dt to 0.6 (default) in D_vel0
;2023-08-30 JL  - copy-pasted a block of code for calculating vel and d_vel from stackplot_cursor.pro
;                 old version was returning 0s 


if 	keyword_set(soho)	then arcsec2km = 718d else arcsec2km = 725d
 
 
NV		= n_elements(line_x[*])

vel		= dblarr(NV-1)			; velocity
D_vel		= dblarr(NV-1)			; first-order uncertainty
D_vel0		= dblarr(NV-1)			; first-order uncertainty based on ds=0.75

if not(keyword_set(ds)) then ds = 0.75d
if not(keyword_set(dt)) then dt = 6d

;calculate velocities
;Note: x=scoord, y=tcoord. If /transpose is set, we have to do ^(-1)
for j=0, NV-2, 1 do begin

  vel[j]	= (line_x[j+1] -line_x[j]) / (line_y[j+1] -line_y[j]) / arcsec2km				; km s^(-1)
  if not(keyword_set(no_transpose)) then vel[j] = 1d/vel[j]

  ;uncertainty. Take the position uncertainty to be half the AIA resolution, i.e., 0.75"
  ;             Take the time uncertainty to be half the cadence, i.e., 6 s
  if    (keyword_set(no_transpose)) then begin
;     D_vel[j]	= 1/2d /arcsec2km  *( (line_x[j+1] -line_x[j] +2*ds)/ (line_y[j+1] -line_y[j] -2*dt) $
; 				     -(line_x[j+1] -line_x[j] -2*ds)/ (line_y[j+1] -line_y[j] +2*dt) )
;     D_vel0[j]	= 1/2d /arcsec2km  *( (line_x[j+1] -line_x[j] +1.5)/ (line_y[j+1] -line_y[j] -2*dt) $
; 				     -(line_x[j+1] -line_x[j] -1.5)/ (line_y[j+1] -line_y[j] +2*dt) )
    D_vel[j]	= sqrt(  arcsec2km^2d *2d* ds^2d *(1d/(line_y[j+1]-line_y[j]))^2d 			$
			+2d*dt^2d *((line_x[j+1]-line_x[j])*arcsec2km/(line_y[j+1]-line_y[j])^2d)^2d )
    D_vel0[j]	= sqrt(  arcsec2km^2d *2d* (0.75d)^2d *(1d/(line_y[j+1]-line_y[j]))^2d 			$
			+2d*6d^2d *((line_x[j+1]-line_x[j])*arcsec2km/(line_y[j+1]-line_y[j])^2d)^2d )
 endif else begin
;     D_vel[j]	= 1/2d *arcsec2km  *( (line_y[j+1] -line_y[j] +2*ds)/ (line_x[j+1] -line_x[j] -2*dt) $
; 				     -(line_y[j+1] -line_y[j] -2*ds)/ (line_x[j+1] -line_x[j] +2*dt) )
;     D_vel0[j]	= 1/2d *arcsec2km  *( (line_y[j+1] -line_y[j] +1.5)/ (line_x[j+1] -line_x[j] -2*dt) $
; 				     -(line_y[j+1] -line_y[j] -1.5)/ (line_x[j+1] -line_x[j] +2*dt) )
    D_vel[j]	= sqrt(  arcsec2km^2d *2d* ds^2d *(1d/(line_x[j+1]-line_x[j]))^2d 			$
			+2d*dt^2d *((line_y[j+1]-line_y[j])*arcsec2km/(line_x[j+1]-line_x[j])^2d)^2d )
    D_vel0[j]	= sqrt(  arcsec2km^2d *2d* (0.75d)^2d *(1d/(line_x[j+1]-line_x[j]))^2d 			$
			+2d*6d^2d *((line_y[j+1]-line_y[j])*arcsec2km/(line_x[j+1]-line_x[j])^2d)^2d )
  endelse

  if not(keyword_set(quiet)) 	then print, '; Velocity for segment ', trim(j+1), ' is ',		$
					       string(vel[j], format='(F9.2)'), ' km/s   +- ',		$
					       string(D_vel[j], format='(F9.2)'), ' km/s    (', 	$
					       string(D_vel0[j], format='(F9.2)'), ' km/s)'
endfor


END