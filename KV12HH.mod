NEURON {
	SUFFIX  KV12
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
	celsius      (degC)
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
	UNITSOFF
	minf= 1/(1+exp(-(v+ 21)/11.3943))
	taum= 150/(1+exp((v+67.56)/34.1479))
	hinf= 1/(1+exp((v+ 22)/11.3943))
	tauh= 15000/(1+exp(-(v+46.56)/44.1479))
	UNITSON
}