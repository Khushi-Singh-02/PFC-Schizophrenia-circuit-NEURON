NEURON {
	SUFFIX  Ca33T
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
	gbar = 0.00001 (S/cm2)
	v 	(mV)
	cai0 =1.84e-4 (mM)
	cao0 = 1.5 (mM)
	celsius =37 (degC)
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
	minf
	taum (ms)
	hinf
	tauh (ms)
	cao (mM)
	cai (mM)
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
	g= gbar* m* h
	eca = 0.043125* (celsius + 273.15)* log(cao / cai)
	curr= 1e-6*g* (v-eca) 			
	ica=curr
}

DERIVATIVE states {
	m'= (minf-m)/taum
	h'= (hinf-h)/tauh
}

PROCEDURE rates(v (mV)){
	UNITSOFF
	minf= 1/(1+exp(-(v+ 45.454426)/5.073))
	taum= 3.394938 + (54.187616/(1+exp((v+40.04)/4.1104)))
	hinf= 1/(1+exp((v+ 74.031965)/8.416382))
	tauh= 109.7 + (0.003816*exp(-v/4.781719))
	UNITSON
}