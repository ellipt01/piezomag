#/bin/sh

gmtset ANNOT_FONT_PRIMARY "Helvetica"
gmtset ANNOT_FONT_SIZE_PRIMARY "10p"
gmtset ANNOT_FONT_SECONDARY "Helvetica"
gmtset ANNOT_FONT_SIZE_SECONDARY "8p"
gmtset HEADER_FONT "Helvetica"
gmtset HEADER_FONT_SIZE "10p"
gmtset LABEL_FONT "Helvetica"
gmtset LABEL_FONT_SIZE "10p"

size=9/9

makecpt -T-5/5/0.25 -Z >| mg.cpt

for i in `seq 0 3`; do

	case $i in
	
		0)
		shiftx=1
		shifty=18
		title="X component"
		;;

		1)
		shiftx=10
		shifty=0
		title="Y component"
		;;

		2)
		shiftx=-10
		shifty=-12
		title="Z component"
		;;

		3)
		shiftx=10
		shifty=0
		title="Total force"
		;;
	
	esac

	sed s'/VAL/'$i'/'g params.dat >| infile
	../main/piez -f infile >| res
	xyz2grd -Gres.grd -I0.1/0.1 -R-10/10/-10/10 res

	if [ $i == 0 ]; then
		grdimage res.grd -JX"$size" -R -P -Cmg.cpt -B5nSWe -X"$shiftx" -Y"$shifty" -K >| results.eps
	else
		grdimage res.grd -JX -R -P -Cmg.cpt -B -X"$shiftx" -Y"$shifty" -K -O >> results.eps
	fi

	grdcontour res.grd -JX -R -P -C0.25 -A0.5f10a0 -L-100/-0.25 -Wta -K -O >> results.eps
	grdcontour res.grd -JX -R -P -C0.25 -A0.5f10a0 -L0/100 -K -O >> results.eps

	# title
	echo "0 11.5 14 0 4 TC $title" | pstext -JX -R -N -K -O >> results.eps

done

psscale -Cmg.cpt -D0/-1.5/5/0.5h -B2.5:"(nT)": -O >> results.eps

evince results.eps &

rm -f res res.grd infile mg.cpt
