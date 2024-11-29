NEURON {
  POINT_PROCESS a3b4
  USEION na WRITE ina
  USEION ca WRITE ica
  USEION k WRITE ik
  POINTER C
  RANGE C2, C1, O, D1, D2, D3, D4, D5, curr
  RANGE g, gmax, del, tdur, ica, ina, ik		:rb
  GLOBAL Erev
  GLOBAL rb, Ru, b, a, Rd1, Rd2, Rd3, Rd4, Rd5, Rr1, Rr2, Rr3, Rr4, Rr5, Rb
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
  Erev	= 3.5    (mV)	: reversal potential
  gmax	= 50  (pS)		: maximal conductance
  tdur = 5 (ms) 		: transmitter duration
  del =  0 (ms)

: Rates
  rb = 1.1e2 (/mM /ms) : ACh binding rate
  Ru = 18 (/ms) : ACh unbinding rate
  b= 6.631 (/ms /mM)
  a= 1.584 (/ms)
  Rd1= 0.595 (/ms)
  Rr1= 0.061 (/ms)
  Rd2= 0.073 (/ms)
  Rr2= 0.025 (/ms)
  Rd3= 0.015 (/ms)
  Rr3= 0.038 (/ms)
  Rd4= 0.039 (/ms)
  Rr4= 7e-4 (/ms)
  Rd5= 5.6e-3 (/ms)
  Rr5= 3e-5 (/ms)
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  curr		(nA)		: current = g*(v - Erev)
  g 		(pS)		: conductance
  C (mM)
  Rb (/ms)
  ica (nA)
  ina (nA)
  ik (nA)
}

STATE {
  C2		: unbound
  C1		: ACh bound
  O		: Open state 
  D1		: desensitized state
  D2		: desensitized state
  D3		: desensitized state
  D4		: desensitized state
  D5		: desensitized state
}

INITIAL {
  C2 = 1
  C1 = 0
  D1 = 0
  D2 = 0
  D3 = 0
  D4 = 0
  D5 = 0
  O = 0
}

BREAKPOINT {
  SOLVE kstates METHOD sparse

  g = gmax * O
  curr = (1e-6) * g * (v - Erev)
  ica = 0.819 * curr
  ina = 0.164 * curr
  ik = 0.0164 * curr
}

KINETIC kstates {
   Rb = C * rb  
  ~ C2 <-> C1 	(Rb, Ru)
  ~ C1 <-> O 	(a, b)
  ~ O <-> D1	(Rd1, Rr1)
  ~ O <-> D2	(Rd2, Rr2)
  ~ O <-> D3	(Rd3, Rr3)
  ~ O <-> D4	(Rd4, Rr4)
  ~ O <-> D5	(Rd5, Rr5)
 
  CONSERVE C2+C1+D1+D2+D3+D4+D5+O = 1
 }
