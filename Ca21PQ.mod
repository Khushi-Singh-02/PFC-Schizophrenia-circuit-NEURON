NEURON {
	SUFFIX  Ca21PQ
	USEION ca READ eca WRITE ica
	RANGE ica, gbar, g, curr
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
	gbar = 0.0001 (S/cm2)
	v 	(mV)
	cai0 = 1.84e-4 (mM)
	cao0 =1.5 (mM)
	celsius = 37 (degC)
}

STATE {m}

ASSIGNED {
	eca (mV)
	ica (mA/cm2)
	curr (mA/cm2)
	g (S/cm2)
	cao (mM)
	cai (mM)
}

INITIAL {
	m= alpha(v)/ (alpha(v)+ beta(v))
	cao = cao0
	cai = cai0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m		
	eca = 0.043125* (celsius + 273.15)* log(cao / cai)
	curr= g* (v-eca) 		:1e-6 in your code, why?	
	ica=curr
}

DERIVATIVE states {
	m'= (1-m)* alpha(v) -m*beta(v)
}

FUNCTION alpha (v (mV)) (/ms) {
	UNITSOFF	
	alpha= 8.5/(1+exp(-(v-8)/12.5))
	UNITSON
}

FUNCTION beta (v (mV)) (/ms) {
	UNITSOFF
	beta= 35/(1+exp((74+v)/14.5))
	UNITSON
}