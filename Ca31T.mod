NEURON {
	SUFFIX  Ca31T
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
	cai0 = 1.84e-4 (mM)
	cao0 = 1.5 (mM)
	celsius = 37 (degC)
}

STATE {
	m(1) 
	h(1)
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
	curr= g* (v-eca) 			
	ica=curr
}

DERIVATIVE states {
	m'= (minf-m)/taum
	h'= (hinf-h)/tauh
}

PROCEDURE rates(v (mV)){
	UNITSOFF
	minf= 1/(1+exp(-(v+ 42.921)/5.1632))
	if (v< -10) {
		taum= -0.855809 + (1.493527*exp(-v/27.4142))
	}
	else {
		taum= 1
	}

	hinf= 1/(1+exp((v+ 72.9)/4.575763))
	tauh= 9.987873 + (0.002883*exp(-v/5.598574))
	UNITSON
}