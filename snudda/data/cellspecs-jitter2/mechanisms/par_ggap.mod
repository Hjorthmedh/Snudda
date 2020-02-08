: ggap.mod
: This is a conductance based gap junction model rather
: than resistance because Traub occasionally likes to 
: set g=0 which of course is infinite resistance.
NEURON {
	POINT_PROCESS gGapPar
	RANGE g, i, vgap
	ELECTRODE_CURRENT i
}
PARAMETER { g = 1e-10 (1/megohm) }
ASSIGNED {
	v (millivolt)
	vgap (millivolt)
	i (nanoamp)
}
BREAKPOINT { i = (vgap - v)*g }






