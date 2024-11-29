NEURON {
	POINT_PROCESS GABAB
	POINTER C
	RANGE R, G, g, gmax, lastrelease, tdur, del
	NONSPECIFIC_CURRENT curr
	GLOBAL Cmax, Cdur
	GLOBAL K1, K2, K3, K4, KD, Erev
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	tdur= 5(ms)
	del=0 (ms)
	Cmax	= 0.5	(mM)		: max transmitter concentration
	Cdur	= 0.3	(ms)		: transmitter duration (rising phase)

:
:	From Kfit with long pulse (5ms 0.5mM)
:
	K1	= 0.52	(/ms mM)	: forward binding rate to receptor
	K2	= 0.0013 (/ms)		: backward (unbinding) rate of receptor
	K3	= 0.098 (/ms)		: rate of G-protein production
	K4	= 0.033 (/ms)		: rate of G-protein decay
	KD	= 100			: dissociation constant of K+ channel
	n	= 4			: nb of binding sites of G-protein on K+
	Erev	= -95	(mV)		: reversal potential (E_K)
	gmax		(umho)		: maximum conductance
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	curr 		(nA)		: current = g*(v - Erev)
	g 		(umho)		: conductance
	C		(mM)		: transmitter concentration
	Gn
}


STATE {
	R				: fraction of activated receptor
	G				: fraction of activated G-protein
}


INITIAL {
	C = 0
	R = 0
	G = 0
}

BREAKPOINT {
	SOLVE bindkin METHOD derivimplicit
	Gn = pow(G,n)
	g = gmax * Gn / (Gn+KD)
	curr = g*(v - Erev)
}


DERIVATIVE bindkin {	
	R' = K1 * C * (1-R) - K2 * R
	G' = K3 * R - K4 * G

}