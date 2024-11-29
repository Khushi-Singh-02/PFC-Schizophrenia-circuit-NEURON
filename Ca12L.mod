NEURON {
	SUFFIX  Ca12L
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
	gbar = 0.05 (S/cm2)	:0.00042
	v 	(mV)
	cai0 = 1.84e-4 (mM)
	cao0 = 1.5 (mM)
	celsius = 37 (degC)
	taum= 10 	(ms)	
	tauh= 59 	(ms)
}

STATE {
	m (1) 
	h (1)
}

ASSIGNED {
	eca (mV)
	ica (mA/cm2)
	curr (mA/cm2)
	g (S/cm2)
	cao (mM)
	cai (mM)	
	minf
	hinf
}

INITIAL {
	rates (v)
	m= minf
	h = hinf
	cao = cao0
	cai = cai0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m* m* h
	eca = 0.043125* (celsius + 273.15)* log(cao / cai)
	curr= g* (v-eca) 			
	ica=curr
}

DERIVATIVE states {
	m'= (minf-m)/taum
	h'= (hinf-h)/tauh
}

PROCEDURE rates(v (mV)){
	UNITSOFF
	minf= 1/(1+exp(-(v+ 30)/6))
	hinf= 1/(1+exp((v+ 80)/6.4))
	UNITSON
}