NEURON {
	SUFFIX  KV14
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
	gbar = 0.0015 (S/cm2)
	v 	(mV)
	ki0 = 150 	(mM)
	ko0 = 5	(mM)
	celsius =37 (degC)
	taum= 3    (ms)
	tauh = 119 (ms)
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
	hinf
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
	g= gbar* m* h
	ek = 0.08625* (celsius + 273.15)* log(ko / ki)
	curr= g* (v-ek) 			
	ik=curr
}

DERIVATIVE states {
	m'= (minf-m)/taum
	h'= (hinf-h)/tauh
}

PROCEDURE rates(v (mV)){
	minf= 1/(1+exp(-(v+ 21.7)/16.9))
	hinf= 1/(1+exp((v+ 73.6)/12.8))
}