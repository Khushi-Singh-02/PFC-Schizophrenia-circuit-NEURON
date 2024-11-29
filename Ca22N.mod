NEURON {
	SUFFIX  Ca22N
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
	celsius = 37 (degC)
	cai0 = 1.84e-4 (mM)
	cao0 = 1.5(mM)
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
	cao (mM)
	cai (mM)
}

INITIAL {
	m= alpham(v)/ (alpham(v)+ betam(v))
	h= alphah(v)/ (alphah(v)+ betah(v))
	cao = cao0
	cai = cai0
}

BREAKPOINT {
	SOLVE states METHOD cnexp		
	g= gbar* m* m* h
	eca = 0.043125* (celsius + 273.15) * log(cao / cai)  
	curr= g* (v-eca) 			
	ica=curr
}

DERIVATIVE states {
	m'= (1-m)* alpham(v) -m*betam(v)
	h'= (1-h)* alphah(v) -h*betah(v)
}

FUNCTION alpham (v (mV)) (/ms) {
	UNITSOFF
	if (v-20> 1e-6) {	
		alpham= 0.1*(v-20)/(1-exp(-(v-20)/10))
	}
	else {
		alpham= 1/(exp(-(v-20)/10))
	}
	UNITSON
}

FUNCTION betam (v (mV)) (/ms) {
	UNITSOFF
	betam= 0.4*exp(-(25+v)/18)
	UNITSON
}

FUNCTION alphah (v (mV)) (/ms) {
	UNITSOFF	
	alphah= 0.01*exp(-(v+50)/10)
	UNITSON
}

FUNCTION betah (v (mV)) (/ms) {
	UNITSOFF
	betah= 0.1/(1+exp(-(17+v)/17))
	UNITSON
}