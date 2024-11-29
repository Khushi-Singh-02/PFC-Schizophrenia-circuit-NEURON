NEURON {
	SUFFIX  KV31
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
	gbar = 0.002 (S/cm2)
	v 	(mV)
	ki0 = 150 	(mM)
	ko0 = 5	(mM)
	celsius = 37 (degC) 
}

STATE {
	m (1)
	h (1)
}

ASSIGNED {
	ek (mV)
	ik (mA/cm2)
	curr (mA/cm2)
	g (S/cm2)
	ki (mM)
	ko (mM)
}

INITIAL {
	m= alpham(v)/ (alpham(v)+ betam(v))
	h= alphah(v)/ (alphah(v)+ betah(v))
	ki = ki0
	ko = ko0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m* m* m* m* h
	ek = 0.08625* (celsius + 273.15)* log(ko / ki)
	curr= g* (v-ek) 			
	ik=curr
}

DERIVATIVE states {
	m'= (1-m)* alpham(v) -m*betam(v)
	h'= (1-h)* alphah(v) -h*betah(v)
}

FUNCTION alpham (v (mV)) (/ms) {
	UNITSOFF	
	alpham= 0.11325* exp(0.025*v)
	UNITSON
}

FUNCTION betam (v (mV)) (/ms) {
	UNITSOFF
	betam= 0.00927* exp(-0.01511*v)
	UNITSON
}

FUNCTION alphah (v (mV)) (/ms) {
	UNITSOFF	
	alphah= 0.00173* exp(-0.1942*v)
	UNITSON
}

FUNCTION betah (v (mV)) (/ms) {
	UNITSOFF
	betah= 0.0935* exp(0.0058*v)
	UNITSON
}