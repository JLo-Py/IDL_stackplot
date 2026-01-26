; JL, Jan 25, 2026

; 0. Compile necessary routines and initiate plotting parameters

; .r /add_your_path/IDL_stackplot/plot_cut_spline.pro
; .r /add_your_path/IDL_stackplot/plot_st_stackplot.pro
; .r /add_your_path/IDL_stackplot/st_stackplot.pro
; .r /add_your_path/IDL_stackplot/nearest_map.pro
; .r /add_your_path/IDL_stackplot/stackplot_cursor.pro
; .r /add_your_path/IDL_stackplot/map_cut_cursor.pro

set_plot, 'X'
device, decomposed = 0 
!p.color = 255
!p.background = 255
xsz=8
ysz=8

; 1. Load SJI file (saved as .sav using: read_iris_l2, sjifile, index, data followed by index2map, index, data, sji1330)

pathsji = '/add_your_path/'
file = 'iris_l2_20240307_210913_3644103603_SJI_1400_t000_crop_DROT'
restore, filename= pathsji + file + '.sav', /verbose
sjimap = sji1400p

sjictable = 'SJI_1400'

; 2. Plot map of interest

t0 = '07-Mar-2024 22:10:47'
t0i=nearest_map(sjimap, ref_time=t0)

wdef, 0, xsz*1e2, ysz*1e2
loadct, 0
plotmap, sjimap[t0i], title = sjimap[t0i].time, /log

; 3. Place cut on the map (left mouse click to insert the points, right click to finish) 

map_cut_cursor, cut_x, cut_y

; 4. Produce stackplot, define the time range first. For curved cut use the /spline keyword.

min_time='22:08:00'
max_time='22:15:00'

outxt=st_stackplot(sjimap, cut_x=cut_x, cut_y=cut_y, min_time=min_time, max_time=max_time, /spline)

; 5. Plot the cut atop of the SJI image. For straight cuts use the plot_cut routine. There are many keywords that handle the appearance of the cut.

loadct, 13
plot_cut_spline, outxt, thick=thick, charsize=1.25, charthick=charthick, color=255, force_first=force_first
; plot_cut, outxt, thick=thick, charsize=1.25, charthick=charthick, color=255

; 6. Plot the stackplot (adjust dmin, dmax, and min and max plotting time if needed.)

min_time_xt='22:09:00'
max_time_xt='22:13:00'

set_plot, 'X'

dmin = 7e0
dmax = 7e2

wdef, 1, xsz*150, ysz*100
IRIS_LCT, sjictable

plot_st_stackplot, outxt, dmin=dmin, dmax=dmax, /log, min_time=min_time_xt, max_time=max_time_xt, /xs, xmargin=[10,2.5], ymargin=[4,1], yminor=yminor, thick=thick, charthick=charthick, xtickinterval=60, charsize=charsize, xticklen=0.03, yticklen=0.02, xminor=6, xtitle='Position along cut [arc sec]'

; 7. Measure velocities using linear fits, use 5 lines by default. Adjust spatial/temporal resolution if needed. Same mouse logic as before.

n_l = 5
ds_iris = 0.33
dt_iris = 2

lines_x=dblarr(n_l,2)
lines_y=dblarr(n_l,2)

vels=dblarr(n_l)
d_vels=dblarr(n_l)

for i=0, n_elements(vels)-1 do begin        & $
    delvarx, line_x, line_y                 & $
    stackplot_cursor, stackplot = outxt, min_time = min_time_xt, max_time = max_time_xt,line_x=line_x, line_y=line_y, vel=vel, d_vel=d_vel & $
    lines_x(i,*)=line_x & $
    lines_y(i,*)=line_y & $
    xyouts, line_x(1)-20, line_y(1)+1, ''+string(i)+'', color=0, charsize=2, charthick=2 & $
    vels(i)=vel & $
    d_vels(i)=d_vel  & $
endfor

; 8. Get mean velocity (optional), reformat the output variables, and save the file

mvel=mean(abs(vels))
mvel=string(mvel, format='(f0.1)')

md_vel=1./n_l*sqrt(total(d_vel^2))
md_vel=string(md_vel, format='(f0.1)')

vels = string(vels, format='(f0.1)')
d_vels=abs(d_vels)
d_vels = string(d_vels, format='(f0.1)')

save, filename = '/add_your_path/example_stackplot.sav', cut_x, cut_y, min_time, max_time, min_time_xt, max_time_xt, outxt, dmin, dmax, vels, d_vels, lines_x, lines_y, n_l, mvel, md_vel

