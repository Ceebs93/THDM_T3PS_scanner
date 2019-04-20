set term postscript enhanced color eps font ',18'

set pm3d map corners2color c4 clip1in

filename ='mhmax_HBwithLHClikelihood.dat'

set xrange [90.:1000.]
set yrange [1.:60.]
set dgrid3d 60,183,1

set zrange [0:*]
set cbrange [0:*]

# Find minimum
set output 'tmp.eps'
splot filename u ($4):($5):(($10) > -0.00001 ? ($10) : 1/0) notit w pm3d
min_z = GPVAL_DATA_Z_MIN
plot filename u ($4):(($10) > -0.00001 ? (($10) < min_z+0.0000001 ? ($10) : 1/0): 1/0) notit w p
min_pos_x = GPVAL_DATA_X_MIN
plot filename u ($5):(($10) > -0.00001 ? (($10) < min_z+0.0000001 ? ($10) : 1/0): 1/0) notit w p
min_pos_y = GPVAL_DATA_X_MIN

print min_pos_x, min_pos_y, min_z

set contour
unset surface
set cntrparam bspline

set table '2sigmacontour.dat'
set cntrparam levels discrete 5.99
splot filename u ($4):($5):(($10)) notit w pm3d
unset table

set table '95CL_all_analyses.dat'
set cntrparam levels discrete 1.0
splot filename u ($4):($5):($13) notit w pm3d
unset table

reset

filename_out ='mhmax_HBwithLHClikelihood.eps'

set pm3d map corners2color c4 clip1in

set xrange [90.:1000.]
set yrange [1.:60.]

set dgrid3d 60,183,1


set size 0.7,1

set palette rgbformulae 30,31,32

set zrange [0.:20.]
set cbrange [0.:20.]

set grid

set xlabel 'M_A [GeV]'
set ylabel 'tan{/Symbol b}'

set output filename_out

set multiplot

splot filename u ($4):($5):(($10)) notit w pm3d

set size 0.545,0.724
set origin 0.081,0.156

unset xtics
unset ytics
unset xlabel
unset ylabel
unset clabel
unset label 1
unset label 2
unset label 3
unset label 4
unset surface
set cntrparam bspline
unset colorbox

set key Left reverse box width -13.4 at 655, 60.0 opaque font ',12'

set label 1 'q@_{MSSM}^{obs}' at 1020,65
set label 2 'm_h^{max} scenario' at 120,64

plot 'mhmax_CMS_obs_wo_SMHiggs_14029.dat' u ($1):($2) w l lt 2 lw 6 lc rgb '#ADFF2F' title '95% CL excl. (CMS)',\
     '2sigmacontour.dat' u ($1):($2) w l lt 1 lw 6 lc rgb 'orange' title 'q@_{MSSM}^{obs}     = 5.99 (reconstr.)'
	 
unset multiplot	 