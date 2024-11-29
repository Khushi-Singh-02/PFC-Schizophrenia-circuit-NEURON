NEURON {
	SUFFIX  KV72
	USEION k READ ek WRITE ik		
	RANGE ik, gbar, g, curr
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(uS) = (microsiemens)
	PI = (pi) (1)
	F = (faraday) (coulomb)
	R = (mole k)  (mV-coulomb/degC)
}

PARAMETER {
	gbar = 0.002(S/cm2)
	v 	(mV)
	ki0 = 150 	(mM)
	ko0 = 5	(mM)
	celsius = 37 (degC) 
}

STATE {m}

ASSIGNED {
	ek (mV)
	ik (mA/cm2)
	curr (mA/cm2)
	g (S/cm2)
	ki (mM)
	ko (mM)
}

INITIAL {
	m= alpha(v)/ (alpha(v)+ beta(v))
	ki = ki0
	ko = ko0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m
	ek = 0.08625* (celsius + 273.15)* log(ko / ki)
	curr= g* (v-ek) 			
	ik=curr
}

DERIVATIVE states {
	m'= (1-m)* alpha(v) -m*beta(v)
}

FUNCTION alpha (v (mV)) (/ms) {
	UNITSOFF	
	alpha= 0.02/(1+exp((40-v)/5))
	UNITSON
}

FUNCTION beta (v (mV)) (/ms) {
	UNITSOFF
	beta= 0.01*exp((17-v)/18)
	UNITSON
}