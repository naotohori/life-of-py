proc vmd_draw_arrow {mol start delta} {
	set end [vecadd $start [vecscale 100 $delta]]
	set middle [vecadd $start [vecscale 0.95 [vecsub $end $start]]]
      graphics $mol color red
	graphics $mol cylinder $start $middle radius 1.0
	graphics $mol cone $middle $end radius 1.50
}

set sel [atomselect top "(resid 99 or resid 144) and name P"]
set xyz [$sel get {x y z}]
set v { {0.90469763 0.4216146  0.06134595} {0.73227063  -0.66028476 0.16674459} }

set vl_max 0
foreach b $v {
	set vl [veclength $b]
	if { $vl_max <= $vl } {
		set vl_max $vl
	}
}

#set vl_min [expr $vl_max / 10.0]
set vl_min 0

set dr_xyz {}
set dr_v {}
foreach a $xyz b $v {
	set vl [veclength $b]
	if { $vl >= $vl_min } {
		lappend dr_xyz $a
		lappend dr_v $b
	}
}

foreach a $dr_xyz b $dr_v {
	draw arrow $a $b
}
