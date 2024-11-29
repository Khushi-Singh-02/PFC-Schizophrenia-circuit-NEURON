NEURON {
    POINT_PROCESS AMPA 
    RANGE R, gmax, g, ina, ik
    USEION na WRITE ina
    USEION k WRITE ik
    POINTER C 
    RANGE g, gmax, rb1, rb2, rd1d2, del, tdur, nai0, nao0
    RANGE C3, C2, C1, D1,D2, D3, O, curr
    GLOBAL Erev
    GLOBAL Rb1, Ru1, Rb2, Ru2, a, b, Rd1, Rr1, Rd2, Rr2, Rd3, Rr3, Rd1d2, Rd2d1, Rd2d3, Rd3d2
} 

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (umho) = (micromho)
    (mM) = (milli/litre)
    (uM) = (micro/liter)
    (mA) = (milliamp)
}

PARAMETER {
    Erev = 0 (mV)	: reversal potential
    gmax = 500 (pS)	: maximal conductance 
    tdur = 5 (ms)
    del = 0  (ms)	:delay
    :nai0 = 15 (mM)
    :nao0 = 140 (mM)

:Rates 

  rb1 = 1.8   (/mM /ms)       : binding of glu to C3  			
  Ru1 = 2.4  (/ms)	           : unbinding of glu from C3
  rb2 = 1   (/mM /ms)          : binding of glu to C2
  Ru2 = 10  (/ms)	           : unbinding  of gly from C2
  b   = 5 (/ms)	           : opening 		
  a   = 16   (/ms)	           : closing 
  Rd1 = 0.7 (/ms)               : desensitisation C3-D1
  Rr1 = 0.15 (/ms)                 : resensitization D1-C3
  Rd2 = 0.1 (/ms)              : desensitisation C2-D2
  Rr2 = 0.0021 (/ms)               : resensitization D2-C2 
  Rd3 = 0.3 (/ms)               : desensitisation O-D3
  Rr3 = 0.015 (/ms)                 : resensitization D3-O
  rd1d2= 10 (/mM /ms)	:d1 to d2
  Rd2d1= 1 (/ms)		:d2 to d1
  Rd2d3= 16 (/ms)		:d2 to d3
  Rd3d2= 12 (/ms)		:d3 to d2
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  curr		(nA)		
  g 		(pS)					
  C (mM) 				
  Rb1  (/ms)
  Rb2 (/ms)
  Rd1d2 (/ms)
  ina (nA)
  ik(nA)
}

STATE {		
  C3		: unbound
  C2		: 1 ligand binds
  C1		: 2 ligands bind
  O		: Open state
  D1		: single ligand, desensitized, from C2
  D2       		: desensitized, from C1
  D3		: desensitized, from C1
 }

INITIAL {		
  C3 = 1
  C2 = 0
  C1 = 0
  D1 = 0
  D2 = 0
  D3=0
  O = 0
}

BREAKPOINT {			
  SOLVE kstates METHOD sparse
  g = gmax * O
  curr = 1e-6* g * (v - Erev)
  ina = 0.87*curr
  ik= 0.13*curr
}

KINETIC kstates {			
    Rb1 = rb1 * C
    Rb2 = rb2 * C
    Rd1d2= rd1d2 * C

  ~ C3 <-> C2	(Rb1, Ru1)			
  ~ C2 <-> C1	(Rb2, Ru2)
  ~ C1 <-> O	(b, a)
  ~ C2 <-> D1	(Rd1, Rr1)
  ~ C1 <-> D2    	(Rd2, Rr2)
  ~ O <-> D3    	(Rd3, Rr3)
  ~ D1 <-> D2    	(Rd1d2, Rd2d1)
  ~ D2 <-> D3    	(Rd2d3, Rd3d2)
 
  CONSERVE C3+C2+C1+D1+D2+D3+O = 1
 }


