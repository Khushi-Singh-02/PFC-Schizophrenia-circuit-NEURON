NEURON {
    POINT_PROCESS GABAA  			
    RANGE R, gmax, g
    USEION cl READ ecl WRITE icl VALENCE -1 	
    NONSPECIFIC_CURRENT ihco3
    POINTER C 
    RANGE rb1, rb2, rd1d2, del, tdur, hco3i, hco3o 
    RANGE C3, C2, C1, D1, D2, O1, O2, curr 
    GLOBAL Erev
    GLOBAL Rb1, Ru1, Rb2, Ru2, a1, b1, a2, b2, Rd1, Rr1, Rd2, Rr2, Rd1d2, Rd2d1
} 

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (umho) = (micromho)
    (mM) = (milli/litre)
    (uM) = (micro/liter)
    (mA) = (milliamp)
    F = (faraday) (coulomb)
    R = (mole k)  (mV-coulomb/degC)
}

PARAMETER {
    celsius = 37 (degC)
    gmax = 28 (pS)	: maximal conductance 
    tdur = 5 (ms)
    del = 0 (ms)
    hco3i = 7 (mM)
    hco3o = 20 (mM)
    Pr= 0.3		:permeability of hco3
    dmax= 3.6(/ms)	:max desensitisation rate
    kh= 2 (mM)	:[GABA] at dmax/2
    nh= 2		:hill coefficient
    ecl= -70mV

:Rates 

  rb1 = 30   (/mM /ms)       : binding of GABA to C3			
  Ru1 = 0.09  (/ms)	           : unbinding of GABA from C3
  rb2 = 15   (/mM /ms)          : binding of GABA to C2
  Ru2 = 0.18  (/ms)	           : unbinding  of GABA from C2
  b1   = 0.2 (/ms)	           : opening (1 ligand, from C2) 		
  a1   = 1.1   (/ms)	           : closing 
  b2   = 2.5 (/ms)	           : opening (2 lignads, from C1) 		
  a2   = 0.142  (/ms)	           : closing 
  Rd1 = 0.006 (/ms)               : desensitisation C3-D1
  Rr1 = 3.5e-4 (/ms)                 : resensitization D1-C3
:Rd2 accounted for later
  Rr2= 0.05 (/ms)
  rd1d2= 1e-8 (/mM /ms)	:d1 to d2
  Rd2d1= 1e-3 (/ms)		:d2 to d1
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  Erev		(mV)	
  curr		(nA)		
  g 		(pS)					:[GABA] 
  ehco3 (mV)			
  C (mM) 				:pointer to give transmitter conc
  Rb1  (/ms)
  Rb2 (/ms)
  Rd2 (/ms)
  Rd1d2(/ms)
  icl (nA)
  ihco3 (nA)
}

STATE {		
  C3		: unbound
  C2		: 1 ligand binds
  C1		: 2 ligands bind
  O1		: Open state, 1 ligand
  O2		: Open state, 2 ligands
  D1		: 1 ligand, desensitized, from C2
  D2       		: 2 ligands, desensitized, from C1
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
  ehco3= log(hco3i/hco3o)*(celsius+ 273.15)*R/F
  Erev= ehco3+ecl
  g = gmax * (O1+O2)
  curr = 1e-6* g * (v - Erev)  
  ihco3 = Pr*curr	
  icl = (1-Pr)*curr
}

KINETIC kstates {			
    Rb1 = rb1 * C
    Rb2 = rb2 * C
    Rd2= dmax/(1+ (kh/C)^nh)
    Rd1d2= rd1d2 * C

  ~ C3 <-> C2	(Rb1, Ru1)			
  ~ C2 <-> C1	(Rb2, Ru2)
  ~ C2 <-> O1	(b1, a1)
  ~ C1 <-> O2	(b2, a2)
  ~ C2 <-> D1	(Rd1, Rr1)
  ~ C1 <-> D2    	(Rd2, Rr2)
  ~ D1 <-> D2    	(Rd1d2, Rd2d1)
 
  CONSERVE C3+C2+C1+D1+D2+O1+O2 = 1
 }