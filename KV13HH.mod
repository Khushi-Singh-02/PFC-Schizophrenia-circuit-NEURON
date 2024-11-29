NEURON {
	SUFFIX  KV13
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
	minf= 1/(1+exp(-(v+ 14.1)/10.3))
	if (v<50) {
		taum= (-0.284*v)+19.16
	}
	else {
		taum= 5
	}

	hinf= 1/(1+exp((v+ 33)/3.7))
	if (v<80) {
		tauh= (-13.76*v)+ 1162.4
	}
	else {
		tauh= 60
	}
	UNITSON
}