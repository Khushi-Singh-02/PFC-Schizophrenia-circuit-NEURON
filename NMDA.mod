NEURON {
  POINT_PROCESS NMDA		
  USEION na WRITE ina 
  USEION ca WRITE ica
  USEION k WRITE ik 
  POINTER C
  RANGE C3, C2, C1, D1, D2, O1, O2, curr, ica, ina, ik
  RANGE g, gmax, rb1, rb2, tr2, del, tdur, nai0, nao0, cai0, cao0, Mg, K0, delta
  GLOBAL Erev
  GLOBAL Rb1, Ru1, Ru2, Rb2, a, b, Rr1, Rr2, Rd1, Rd2, Ro1, Ro2 
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (umho) = (micromho)
  (mM) = (milli/liter)
  (uM) = (micro/liter)
  (mA) = (milliamp)
  R= (mole k) (mV-coulomb/degC)
  F= (faraday) (coulomb)
}
:check all values
PARAMETER {
  Erev	= 0    (mV)		: reversal potential
  gmax	=  48 (pS)	: maximal conductance, too high?
  tdur = 5 (ms)
  del =  0 (ms)
  tr2= 0.1  (mM)		:gly conc const
  :nai0 = 15 (mM)
  :nao0 = 140 (mM)
  :cai0 = 1.84e-4 (mM)
  :cao0 = 1.5 (mM)

: Rates

  rb1 = 0.257   (/mM /ms)   : binding of glu to C3  			
  Ru1 = 0.06  (/ms)	           : unbinding of glu from C3
  rb2 = 0.775   (/mM /ms)   : binding of gly to C2
  Ru2 = 1.434  (/ms)	           : unbinding  of gly from C2
  b   = 3.67 (/ms)	           : opening 		
  a   = 0.72   (/ms)	      : closing 
  Ro2 = 1.62(/ms)		: O1-O2
  Ro1 = 0.74(/ms)		: O2-O1
  Rd1 = 0.008 (/ms)          : desensitisation C3-D1
  Rr1 = 5e-4 (/ms)           : resensitization D1-C3
  Rd2 = 0.006 (/ms)          : desensitisation C2-D2
  Rr2 = 0.038 (/ms)           : resensitization D2-C2

:For Mg
Mg=1(mM)	:extracellular Mg conc
K0=4.1(mM)	:IC50 at 0 mV
delta= 0.8		:electrical distance
T= 37 (degC)
z=2		:Mg charge
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  curr		(nA)		
  g 		(pS)
  C (mM) 				:pointer to pass conc of transmitter 
  Rb1  (/ms)
  Rb2 (/ms)
  ica (nA)
  ina (nA)
  ik (nA)
  MB (1)
}

STATE {		
  C3		: unbound, traps glu
  C2		: glu bound, traps gly, controls channel opening
  C1		: gly bound, channel opens
  O1		: Open state 1 
  O2		: Open state 2 (kinetically distinct)
  D1		: desensitized, from C3
  D2       		: desensitized, from C2
 }

INITIAL {		
  C3 = 1
  C2 = 0
  C1 = 0
  D1 = 0
  D2 = 0
  O1 = 0
  O2 = 0
}




BREAKPOINT {			
  SOLVE kstates METHOD sparse
  MgBlock(v) 
  g = gmax * (O1+O2)
  curr = (1e-6) *MB* g * (v - Erev)
  ica = 0.816*curr
  ina = 0.15*curr
  ik = 0.06*curr	
}


KINETIC kstates {			
    Rb1 = rb1 * C
    Rb2 = rb2 * tr2

  ~ C3 <-> C2	(Rb1, Ru1)			
  ~ C2 <-> C1	(Rb2, Ru2)
  ~ C1 <-> O1	(b, a)
  ~ O1 <-> O2	(Ro2, Ro1)
  ~ C3 <-> D1	(Rd1, Rr1)
  ~ C2 <-> D2    	(Rd2, Rr2)
 
  CONSERVE C3+C2+C1+D1+D2+O1+O2 = 1
 }

FUNCTION MgBlock(v(mV)){ 
  MB= 1/(1+(Mg/K0)*exp((-z)*delta*F*v/R/(T+273.15)))
}

