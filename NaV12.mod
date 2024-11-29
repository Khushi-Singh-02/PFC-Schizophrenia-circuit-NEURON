NEURON {
	SUFFIX NaV12 	
	USEION na READ ena WRITE ina
	RANGE gmax, ina, g
	RANGE C1, C2, O1, O2, I1, I2
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	nai0 = 15 (mM)
	nao0 = 140 (mM)
	celsius = 37 (degC)
	gbar  = 0.2 (S/cm2) 
	
	C1C2b2	  = 16
	C1C2v2    = -5
	C1C2k2	  = -10
	
	C2C1b1	  = 3
	C2C1v1    = -35
	C2C1k1	  = 10
	C2C1b2	  = 16
	C2C1v2    = -5
	C2C1k2	  = -10

	C2O1b2	  = 16
	C2O1v2    = -10
	C2O1k2	  = -10
	
	O1C2b1	  = 3
	O1C2v1    = -40
	O1C2k1	  = 10
	O1C2b2	  = 16
	O1C2v2    = -10
	O1C2k2	  = -10
	
	C2O2b2	  = 0.13
	C2O2v2	  = -20
	C2O2k2	  = -15
	
	O2C2b1	  = 2
	O2C2v1	  = -60
	O2C2k1	  = 6
	O2C2b2	  = 0.7
	O2C2v2	  = -10
	O2C2k2	  = -15
	
	O1I1b1	  = 3
	O1I1v1	  = -41
	O1I1k1	  = 12
	O1I1b2	  = 16
	O1I1v2	  = -11
	O1I1k2	  = -12
	
	I1O1b1	  = 1e-5
	I1O1v1	  = -42
	I1O1k1	  = 10
	
	I1C1b1	  = 0.55
	I1C1v1	  = -65
	I1C1k1	  = 10
	
	C1I1b2	  = 0.55
	C1I1v2	  = -65
	C1I1k2	  = -11
	
	I1I2b2	  = 0.0022
	I1I2v2	  = -90
	I1I2k2	  = -5

	I2I1b1	  = 0.017
	I2I1v1	  = -90
	I2I1k1	  = 15
	
}

ASSIGNED {
	ena (mV)
	ina  (mA/cm2)
	g   (mho/cm2)
	nai (mM)
	nao (mM)	
	RC1C2 (/ms)
	RC2C1 (/ms)
	RC2O1 (/ms)
	RO1C2 (/ms)
	RC2O2 (/ms)
	RO2C2 (/ms)
	RO1I1 (/ms)
	RI1O1 (/ms)
	RI1I2 (/ms)
	RI2I1 (/ms)
	RI1C1 (/ms)
	RC1I1 (/ms)
	Q10 (1)
}

STATE {
	C1
	C2
	O1
	O2
	I1
	I2
}


INITIAL {
	Q10 = 3^((celsius-20 (degC))/10 (degC))
	:SOLVE kin
	: STEADYSTATE sparse ???
	C1=1
	C2=0
	O1=0
	O2=0
	I1=0
	I2=0
	nai=nai0
	nao=nao0
}

BREAKPOINT {
	SOLVE kstates METHOD sparse
	g = gbar * (O1 + O2)	
	ena= 0.08625* (celsius + 273.15)* log(nao / nai)
	ina = g * (v - ena)   	
}

KINETIC kstates {
	rates(v)
	
	~ C1 <->  C2 	(RC1C2, RC2C1)
	~ C2 <->  O1 	(RC2O1, RO1C2)
	~ C2 <->  O2 	(RC2O2, RO2C2)
	~ O1 <->  I1 	(RO1I1, RI1O1)
	~ I1   <->  C1 	(RI1C1, RC1I1)
	~ I1   <->  I2 	(RI1I2, RI2I1)
	
	CONSERVE O1 + O2 + C1 + C2 + I1 + I2 = 1
}

FUNCTION rates2(v, b, vv, k) {
	rates2 = (b/(1+exp((v-vv)/k)))
}

PROCEDURE rates(v(mV)) {
UNITSOFF
	RC1C2 = Q10*(rates2(v, C1C2b2, C1C2v2, C1C2k2))
	RC2C1 = Q10*(rates2(v, C2C1b1, C2C1v1, C2C1k1) + rates2(v, C2C1b2, C2C1v2, C2C1k2))
	RC2O1 = Q10*(rates2(v, C2O1b2, C2O1v2, C2O1k2))
	RO1C2 = Q10*(rates2(v, O1C2b1, O1C2v1, O1C2k1) + rates2(v, O1C2b2, O1C2v2, O1C2k2))
	RC2O2 = Q10*(rates2(v, C2O2b2, C2O2v2, C2O2k2))
	RO2C2 = Q10*(rates2(v, O2C2b1, O2C2v1, O2C2k1) + rates2(v, O2C2b2, O2C2v2, O2C2k2))
	RO1I1 = Q10*(rates2(v, O1I1b1, O1I1v1, O1I1k1) + rates2(v, O1I1b2, O1I1v2, O1I1k2))
	RI1O1 = Q10*(rates2(v, I1O1b1, I1O1v1, I1O1k1))
	RI1C1 = Q10*(rates2(v, I1C1b1, I1C1v1, I1C1k1))
	RC1I1 = Q10*(rates2(v, C1I1b2, C1I1v2, C1I1k2))
	RI1I2 = Q10*(rates2(v, I1I2b2, I1I2v2, I1I2k2))
	RI2I1 = Q10*(rates2(v, I2I1b1, I2I1v1, I2I1k1))
UNITSON
}
