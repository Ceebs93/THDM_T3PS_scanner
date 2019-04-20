set term postscript enhanced color font ',18' eps

filename_out = '2Higgses.eps'

filename1 ='pdf1.dat'
filename2 ='pdf2.dat'
filename3 ='pdf3.dat'

# First, do the contour lines

set xrange [122.:128.]
set yrange [122.:128.]

set dgrid3d 31,31,1
set pm3d map corners2color c1 clip4in

# Find minimum
set output 'tmp/tmp.eps'
splot filename1 u ($1):($2):($9) notit w pm3d
min_z_1 = GPVAL_DATA_Z_MIN
plot filename1 u ($1):($9 < min_z_1+0.000001 ? $9 : 1/0) notit w p
min_pos_x_1 = GPVAL_DATA_X_MIN
plot filename1 u ($2):($9 < min_z_1+0.00001 ? $9 : 1/0) notit w p
min_pos_y_1 = GPVAL_DATA_X_MIN

print " "
print "      -------------- "
print "      best-fit point (pdf = 1):"
print "      -------------- "
print "      minimal Chi^2 Value = ", min_z_1
print "      (mh1, mh2) = (", min_pos_x_1,",", min_pos_y_1,")"

# Find minimum
set output 'tmp/tmp.eps'
splot filename2 u ($1):($2):($9) notit w pm3d
min_z_2 = GPVAL_DATA_Z_MIN
plot filename2 u ($1):($9 < min_z_2+0.000001 ? $9 : 1/0) notit w p
min_pos_x_2 = GPVAL_DATA_X_MIN
plot filename2 u ($2):($9 < min_z_2+0.00001 ? $9 : 1/0) notit w p
min_pos_y_2 = GPVAL_DATA_X_MIN

print "      -------------- "
print "      best-fit point (pdf = 2):"
print "      -------------- "
print "      minimal Chi^2 Value = ", min_z_2
print "      (mh1, mh2) = (", min_pos_x_2,",", min_pos_y_2,")"


# Find minimum
set output 'tmp/tmp.eps'
splot filename3 u ($1):($2):($9) notit w pm3d
min_z_3 = GPVAL_DATA_Z_MIN
plot filename3 u ($1):($9 < min_z_3+0.000001 ? $9 : 1/0) notit w p
min_pos_x_3 = GPVAL_DATA_X_MIN
plot filename3 u ($2):($9 < min_z_3+0.00001 ? $9 : 1/0) notit w p
min_pos_y_3 = GPVAL_DATA_X_MIN

print "      -------------- "
print "      best-fit point (pdf = 3):"
print "      -------------- "
print "      minimal Chi^2 Value = ", min_z_3
print "      (mh1, mh2) = (", min_pos_x_3,",", min_pos_y_3,")"
print "      -------------- "


plot filename1 u ($3):($3) notit w p
dm1 = GPVAL_DATA_X_MIN
plot filename1 u ($4):($4) notit w p
dm2 = GPVAL_DATA_X_MIN

plot filename1 u ($5):($5) notit w p
mu1 = GPVAL_DATA_X_MIN
plot filename1 u ($6):($6) notit w p
mu2 = GPVAL_DATA_X_MIN


set contour
unset surface
set cntrparam bspline
unset clip
unset colorbox

set table 'tmp/1sigmacontour_pdf1.dat'
set cntrparam levels discrete 2.2958
splot filename1 u ($1):($2):($9-min_z_1) notit w pm3d
unset table

set table 'tmp/2sigmacontour_pdf1.dat'
set cntrparam levels discrete 5.99
splot filename1 u ($1):($2):($9-min_z_1) notit w pm3d
unset table

set table 'tmp/1sigmacontour_pdf2.dat'
set cntrparam levels discrete 2.2958
splot filename2 u ($1):($2):($9-min_z_2) notit w pm3d
unset table

set table 'tmp/2sigmacontour_pdf2.dat'
set cntrparam levels discrete 5.99
splot filename2 u ($1):($2):($9-min_z_2) notit w pm3d
unset table

set table 'tmp/1sigmacontour_pdf3.dat'
set cntrparam levels discrete 2.2958
splot filename3 u ($1):($2):($9-min_z_3) notit w pm3d
unset table

set table 'tmp/2sigmacontour_pdf3.dat'
set cntrparam levels discrete 5.99
splot filename3 u ($1):($2):($9-min_z_3) notit w pm3d
unset table

reset

set output filename_out
set size 2.2,1.1
set origin 0.0, 0.0

set multiplot

# set size 1.0,1.0
set size 0.7,1
set origin 0.05, 0.1

set pm3d map corners2color c1 clip4in
set palette rgbformulae 30,31,32

set xrange [122.:128.]
set yrange [122.:128.]

set dgrid3d 31,31,1

set xlabel 'm_1 [GeV]' offset 0, -0.6
set ylabel 'm_2 [GeV]' offset -1.2
set xtics 122.,1
set mxtics 10
set ytics 122.,1.
set mytics 10

set zrange [0:20]
set cbrange [0:20]

# set multiplot
set grid

set label 1 point ps 1.5 pt 3 lc rgb 'green' at min_pos_x_1, min_pos_y_1 front
set label 2 '{/Symbol m}_1 = '.gprintf("%4.2f",mu1).', {/Symbol m}_2 = '.gprintf("%4.2f",mu2).', {/Symbol D}m_1 = '.gprintf("%4.2f",dm1).' GeV, {/Symbol D}m_2 = '.gprintf("%4.2f",dm2).' GeV' font ',13' at 122.1,128.2 front
set label 3 '(a) Box pdf.' at 122.0,120.0 font ',22' front
set label 4 '{/Symbol D}{/Symbol c}_{tot}^2' at 128.1,128.5 front

splot filename1 u ($1):($2):($9-min_z_1+0.0001) notit w pm3d

set size 0.545,0.7221
set origin 0.132,0.258

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

set key Left reverse spacing 1.2 box width -2 font ',14' at 128, 128 opaque

plot 'tmp/1sigmacontour_pdf1.dat' u ($1):($2) w l lt 1 lw 4 lc rgb '#DDDDDD' title '68% C.L.',\
	 'tmp/2sigmacontour_pdf1.dat' u ($1):($2) w l lt 2 lw 4 lc rgb '#DDDDDD' title '95% C.L.'

reset
	 
set size 0.7,1
set origin 0.75, 0.1

set pm3d map corners2color c1 clip4in
set palette rgbformulae 30,31,32

set xrange [122.:128.]
set yrange [122.:128.]

set dgrid3d 31,31,1

set xlabel 'm_1 [GeV]' offset 0, -0.6
set ylabel 'm_2 [GeV]' offset -1.2
set xtics 122.,1
set mxtics 10
set ytics 122.,1.
set mytics 10

set zrange [0:20]
set cbrange [0:20]

set grid

set label 21 point ps 1.5 pt 3 lc rgb 'green' at min_pos_x_2, min_pos_y_2 front
set label 22 '{/Symbol m}_1 = '.gprintf("%4.2f",mu1).', {/Symbol m}_2 = '.gprintf("%4.2f",mu2).', {/Symbol D}m_1 = '.gprintf("%4.2f",dm1).' GeV, {/Symbol D}m_2 = '.gprintf("%4.2f",dm2).' GeV' font ',13' at 122.1,128.2 front
set label 23 '(b) Gaussian pdf.' at 122.0,120.0 font ',22' front
set label 24 '{/Symbol D}{/Symbol c}_{tot}^2' at 128.1,128.5 front

splot filename2 u ($1):($2):($9-min_z_2+0.0001) notit w pm3d

set size 0.545,0.7221
set origin 0.832,0.258

unset xtics
unset ytics
unset xlabel
unset ylabel
unset clabel
unset label 21
unset label 22
unset label 23
unset label 24
unset surface
set cntrparam bspline
unset colorbox

set key Left reverse spacing 1.2 box width -2 font ',14' at 128, 128 opaque

plot 'tmp/1sigmacontour_pdf2.dat' u ($1):($2) w l lt 1 lw 4 lc rgb '#DDDDDD' title '68% C.L.',\
	 'tmp/2sigmacontour_pdf2.dat' u ($1):($2) w l lt 2 lw 4 lc rgb '#DDDDDD' title '95% C.L.'	 


reset
	 
set size 0.7,1
set origin 1.45, 0.1

set pm3d map corners2color c1 clip4in
set palette rgbformulae 30,31,32

set xrange [122.:128.]
set yrange [122.:128.]

set dgrid3d 31,31,1

set xlabel 'm_1 [GeV]' offset 0, -0.6
set ylabel 'm_2 [GeV]' offset -1.2
set xtics 122.,1
set mxtics 10
set ytics 122.,1.
set mytics 10

set zrange [0:20]
set cbrange [0:20]

set grid

set label 31 point ps 1.5 pt 3 lc rgb 'green' at min_pos_x_3, min_pos_y_3 front
set label 32 '{/Symbol m}_1 = '.gprintf("%4.2f",mu1).', {/Symbol m}_2 = '.gprintf("%4.2f",mu2).', {/Symbol D}m_1 = '.gprintf("%4.2f",dm1).' GeV, {/Symbol D}m_2 = '.gprintf("%4.2f",dm2).' GeV' font ',13' at 122.1,128.2 front
set label 33 '(c) Box+Gaussian pdf.' at 122.0,120.0 font ',22' front
set label 34 '{/Symbol D}{/Symbol c}_{tot}^2' at 128.1,128.5 front

splot filename3 u ($1):($2):($9-min_z_3+0.0001) notit w pm3d

set size 0.545,0.7221
set origin 1.532,0.258

unset xtics
unset ytics
unset xlabel
unset ylabel
unset clabel
unset label 31
unset label 32
unset label 33
unset label 34
unset surface
set cntrparam bspline
unset colorbox

set key Left reverse spacing 1.2 box width -2 font ',14' at 128, 128 opaque

plot 'tmp/1sigmacontour_pdf3.dat' u ($1):($2) w l lt 1 lw 4 lc rgb '#DDDDDD' title '68% C.L.',\
	 'tmp/2sigmacontour_pdf3.dat' u ($1):($2) w l lt 2 lw 4 lc rgb '#DDDDDD' title '95% C.L.'	 
	 
	 	 
unset multiplot	 

