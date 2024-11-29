NEURON{
	SUFFIX leak
	NONSPECIFIC_CURRENT curr
	RANGE curr, il, el, g
}

UNITS{
	(mA) = (milliamp)
	(mS) = (millimho)
	(mV) = (millivolt)
}

PARAMETER{
	v (mV)
	gl = 0.3 (mS/cm2)
	el = -65 (mV)
}

ASSIGNED{
	curr (mA/cm2)
	il (mA/cm2)
}

BREAKPOINT{
	il= (1e-3)*gl*(v-el)
}