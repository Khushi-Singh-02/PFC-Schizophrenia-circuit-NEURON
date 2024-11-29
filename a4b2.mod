NEURON {
  POINT_PROCESS a4b2
  USEION na WRITE ina
  USEION ca WRITE ica
  POINTER C
  RANGE C0, C1, C2, D1, D2, O, D0, curr, ica, ina
  RANGE g, gmax, rb, rdb1, rdb2, del, tdur
  GLOBAL Erev
  GLOBAL Rb, Ru, Rdb1, Rdb2, Rdu1, Rdu2, a, b, Rrc, Rdc 
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (umho) = (micromho)
  (mM) = (milli/liter)
  (uM) = (micro/liter)
  (mA) = (milliamp)
}

PARAMETER {
  Erev	= 0    (mV)	: reversal potential
  gmax	= 48  (pS)	: maximal conductance
  tdur	= 5 (ms)
  del 	= 0 (ms)

: Rates

  rb = 1.1e2   (/mM /ms)   : binding  			
  Ru = 18  (/ms)	           : unbinding
  b = 50 (/ms)	           : opening 		
  a = 1.2   (/ms)	      : closing 
  Rd = 0.004 (/ms)          : desensitization 
  Rr = 1e-5 (/ms)           : resensitization
  Rdu1 = 0.046(/ms)	       : unbinding from 1st desensitization state
  rdb1 = 5.7e2    (/mM /ms)	: binding to 1st desensitization state
  Rdu2 = 0.023   (/ms)	  : unbinding from 2nd desensitization state
  rdb2 = 1.1e3 (/mM /ms)    : binding to 2nd desensitization state
  Rrc = 0.004 (/ms)          : resensitization to closed state
  Rdc = 5e-6 (/ms)           : desensitization from closed state
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  curr	(nA)			: current = g*(v - Erev)
  g 		(pS)		: conductance
  C (mM)
  Rb		(/ms)    
  Rdb1  (/ms)
  Rdb2 (/ms)
  ica (nA)
  ina (nA)
}

STATE {
  C0		: unbound
  C1		: single ACh bound
  C2		: double ACh bound
  O		: Open state 
  D0		: double ACh bound, desensitized
  D1       		: single ACh bound, desensitized
  D2 		: unbound, desensitized
}

INITIAL {
  C0 = 1
  C1 = 0
  C2 = 0
  D0 = 0
  D1 = 0
  D2 = 0
  O = 0
}

BREAKPOINT {
  SOLVE kstates METHOD sparse

  g = gmax * O
  curr = (1e-6) * g * (v - Erev)
  ica = 0.6*curr
  ina = 0.4*curr
}

KINETIC kstates {
    Rb = rb * C
    Rdb1 = rdb1 * C
    Rdb2 = rdb2 * C

  ~ C0  <-> C1	(2*Rb, Ru)
  ~ C1 <-> C2	(Rb, 2*Ru)
  ~ C2 <-> O	(b, a)
  ~ O <-> D0	(Rd, Rr)
  ~ D0 <-> D1	(Rdu1, Rdb1)
  ~ D1 <-> D2    (Rdu2, Rdb2)
  ~ D2 <-> C0    (Rrc, Rdc)
 
  CONSERVE C0+C1+C2+D0+D1+D2+O = 1
 }