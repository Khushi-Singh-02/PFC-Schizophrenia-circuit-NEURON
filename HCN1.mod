NEURON {
  SUFFIX HCN1
  NONSPECIFIC_CURRENT ih
  RANGE gbar, g, e, v50, htau, hinf
  RANGE gfactor, htaufactor
}
 
UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
}

PARAMETER {
  celsius	(degC)
  gbar = 0.002 (S/cm2)
  e= -30	(mV)
  v50= -73	(mV)
  gfactor = 1
  htaufactor = 4.78
}
 
STATE {
  h
}
 
ASSIGNED {
  ih	  (mA/cm2) 
  hinf
  htau    (ms)
  v	  (mV)
  g       (mho/cm2)
}

PROCEDURE giassign () { 
  g = gbar*h*gfactor
  ih = g*(v-e)
}
 
BREAKPOINT {
  SOLVE states METHOD cnexp
  giassign()
}
 
DERIVATIVE states { 
  rates(v)
  h'= (hinf- h)/ htau
}

INITIAL { 
  rates(v)
  h = hinf
  giassign()
}

PROCEDURE rates(v (mV)) {
  UNITSOFF
  hinf = 1/(1+exp(0.151*(v-v50)))
  htau = htaufactor*exp((0.033*(v+75)))/(0.011*(1+exp(0.083*(v+75))))
  UNITSON
}
