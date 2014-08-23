#/bin/sh

# Name of script:
# calcomp.sh
#
# Require:
# GMT (Generic Mapping Tools)
#
# Description:
# This script calculates seismomagnetic field due to strike-slip, dip-slip or tensile opening fault motion
# using the program "main/piez", and draw contour figures of X, Y, Z component and Total force using GMT.
# The range of calculation is fixed in x(EW)=[-10:0.1:10](km), y(NS)=[-10:0.1:10](km)
# If GMT is not installed on your system, this script will abort with some errors.
#
# Usage:
# $ calcomp.sh <parameter file name>
#

if [ -z "$1" ]; then
	echo "USAGE: $0 <parameter file name>."
	exit
fi

fnparam="$1"
fnout="results.eps"

# GMT settings
gmtset PAGE_ORIENTATION "portrait"
gmtset PAPER_MEDIA "a4"
gmtset ANNOT_FONT_PRIMARY "Helvetica"
gmtset ANNOT_FONT_SIZE_PRIMARY "10p"
gmtset ANNOT_FONT_SECONDARY "Helvetica"
gmtset ANNOT_FONT_SIZE_SECONDARY "8p"
gmtset HEADER_FONT "Helvetica"
gmtset HEADER_FONT_SIZE "10p"
gmtset LABEL_FONT "Helvetica"
gmtset LABEL_FONT_SIZE "10p"

size=9/9

# read fault parameters from input file
L1=`cat example_params.dat | gawk '{if($1 == "flength1"){print $3}}'`
L2=`cat example_params.dat | gawk '{if($1 == "flength2"){print $3}}'`
W1=`cat example_params.dat | gawk '{if($1 == "fwidth1"){print $3}}'`
W2=`cat example_params.dat | gawk '{if($1 == "fwidth2"){print $3}}'`
STRIKE=`cat example_params.dat | gawk '{if($1 == "fstrike"){print $3}}'`
DIP=`cat example_params.dat | gawk '{if($1 == "fdip"){print $3}}'`

U1=`cat example_params.dat | gawk '{if($1 == "u1"){print $3}}'`
U2=`cat example_params.dat | gawk '{if($1 == "u2"){print $3}}'`
U3=`cat example_params.dat | gawk '{if($1 == "u3"){print $3}}'`

LENGTH=`echo $L1+$L2 | bc`
WIDTH=`echo $W1+$W2 | bc`
PI=3.141593
SS=`echo 3.14/2 - $STRIKE*$PI/180 | bc -l`
DD=`echo $DIP*$PI/180 | bc -l`

# create cpt file
makecpt -T-5/5/0.25 -Z >| mg.cpt

# main title
title="FAULT PARAMETERS: dimension (L,W)=($LENGTH, $WIDTH)"
echo "8 48 14 0 4 TC $title" | pstext -JX"$size" -R-10/10/-10/10 -N -K >| $fnout

title="dislocation=($U1, $U2, $U3), strike and dip angle=($STRIKE, $DIP)"
echo "8 46.5 14 0 4 TC $title" | pstext -JX"$size" -R -N -K -O >> $fnout

# draw contours for each components
# X = 0, Y = 1, Z = 2, F = 3
for i in 1 2 3 0; do

	case $i in
	
		1)
		shiftx=-1
		shifty=14
		title="X component"
		;;

		2)
		shiftx=10
		shifty=0
		title="Y component"
		;;

		3)
		shiftx=-10
		shifty=-12
		title="Z component"
		;;

		0)
		shiftx=10
		shifty=0
		title="Total force"
		;;
	
	esac

	sed s'/VAL/'$i'/'g $fnparam >| _tmp_infile_

	# calculate seismomagnetic field in range x=[-10:0.1:10], y=[-10:0.1:10]
	../main/piez -f _tmp_infile_ -r -10/10/-10/10 -i 0.1/0.1 >| res

	# create contour figure
	surface -Gres.grd -I0.1/0.1 -R-10/10/-10/10 -: res

	# draw color contour
	grdimage res.grd -JX -R -P -Cmg.cpt -B5nSWe -X"$shiftx" -Y"$shifty" -K -O >> $fnout

	# draw fault plane
	echo -$L1 -$W2 $SS $DD | gawk '{printf "%f %f\n", $1*cos($3)-$2*cos($4)*sin($3), $1*sin($3)+$2*cos($4)*cos($3)}' >| _tmp_fault_
	echo -$L1  $W1 $SS $DD | gawk '{printf "%f %f\n", $1*cos($3)-$2*cos($4)*sin($3), $1*sin($3)+$2*cos($4)*cos($3)}' >> _tmp_fault_
	echo  $L2  $W1 $SS $DD | gawk '{printf "%f %f\n", $1*cos($3)-$2*cos($4)*sin($3), $1*sin($3)+$2*cos($4)*cos($3)}' >> _tmp_fault_
	echo  $L2 -$W2 $SS $DD | gawk '{printf "%f %f\n", $1*cos($3)-$2*cos($4)*sin($3), $1*sin($3)+$2*cos($4)*cos($3)}' >> _tmp_fault_
	psxy _tmp_fault_ -JX -R -P -W4 -L -K -O >> $fnout

	echo -$L1  $W1 $SS $DD | gawk '{printf "%f %f\n", $1*cos($3)-$2*cos($4)*sin($3), $1*sin($3)+$2*cos($4)*cos($3)}' >| _tmp_fault_
	echo  $L2  $W1 $SS $DD | gawk '{printf "%f %f\n", $1*cos($3)-$2*cos($4)*sin($3), $1*sin($3)+$2*cos($4)*cos($3)}' >> _tmp_fault_
	psxy _tmp_fault_ -JX -R -P -W12 -K -O >> $fnout

	# draw contour lines
	grdcontour res.grd -JX -R -P -C0.25 -A0.5f10a0 -L-100/-0.25 -Wta -K -O >> $fnout
	grdcontour res.grd -JX -R -P -C0.25 -A0.5f10a0 -L0/100 -K -O >> $fnout

	# title
	echo "0 11.5 14 0 4 TC $title" | pstext -JX -R -N -K -O >> $fnout

done

# color scale
psscale -Cmg.cpt -D0/-1.5/5/0.5h -B2.5:"(nT)": -O >> $fnout

evince $fnout &

rm -f res res.grd _tmp_infile_ _tmp_fault_ mg.cpt
