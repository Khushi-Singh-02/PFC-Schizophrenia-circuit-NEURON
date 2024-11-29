NEURON {
	SUFFIX  KV11
	USEION k READ ek WRITE ik		:READ???
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
	gbar = 0.0015(S/cm2)
	v 	(mV)
	ki0 = 150 	(mM)
	ko0 = 5	(mM)
	celsius= 37(degC)
}

STATE {
	m(1) 
	h(1)
}

ASSIGNED {
	ek (mV)
	ik (mA/cm2)
	curr (mA/cm2)
	g (S/cm2)
	minf
	taum (ms)
	hinf
	tauh (ms)
	ki (mM)
	ko (mM)
}

INITIAL {
	rates (v)
	m= minf
	h = hinf
	ki = ki0
	ko = ko0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m* h* h
	ek = 0.08625* (celsius + 273.15)* log(ko / ki)
	curr= g* (v-ek) 			
	ik=curr
}

DERIVATIVE states {
	m'= (minf-m)/taum
	h'= (hinf-h)/tauh
}

PROCEDURE rates(v (mV)){
	UNITSOFF
	minf= 1/(1+exp(-(v+ 30.5)/11.3943))
	taum= 30/(1+exp((v+76.56)/26.1479))
	hinf= 1/(1+exp((v+ 30)/27.3943))
	tauh= 15000/(1+exp(-(v+160.56)/100))
	UNITSON
}