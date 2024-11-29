NEURON {
	SUFFIX  KV32
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
	celsius =37 (degC)
}

STATE {m}

ASSIGNED {
	ek (mV)
	ik (mA/cm2)
	curr (mA/cm2)
	g (S/cm2)
	minf
	taum (ms)
	ki (mM)
	ko (mM)
}

INITIAL {
	rates (v)		
	m= minf
	ki = ki0
	ko = ko0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m* m
	ek = 0.08625* (celsius + 273.15)* log(ko / ki)
	curr= g* (v-ek) 			
	ik=curr
}

DERIVATIVE states {
	m'= (minf-m)/taum
}

PROCEDURE rates(v (mV)){
	minf= 1/(1+exp(-(v+ 0.373267)/8.568187))
	taum= 3.241643+ (19.106496/(1+exp((v- 19.22)/4.451533)))
}