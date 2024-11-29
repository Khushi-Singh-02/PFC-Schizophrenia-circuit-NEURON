NEURON {
  SUFFIX NaKpump
  USEION na READ nai WRITE ina
  USEION k WRITE ik
  RANGE  inapump,ipumpmax,n,km, nai0, nao0 
}

UNITS {
  (molar) = (1/liter)
  (pA) =  (picoamp)
  (mV) =	(millivolt)
  (uS) =  (micromho)
  (mA) =	(milliamp)
  (mM) =	(millimolar)
}

PARAMETER {
  dt  (ms)
  celsius = 37 (degC)
  ipumpmax  = 0.0092 (mA/cm2)
  km = 10.0 (mM)
  n  = 1.5
  nai0 = 5 (mM)
  nao0 = 154 (mM)      
}

ASSIGNED { 
  ina	   (mA/cm2)
  ik	   (mA/cm2)
  nai (mM)
  nao (mM)
  inapump (mA/cm2)
}

INITIAL{
  nai = nai0
  nao = nao0 
}

BREAKPOINT {
  inapump = ipumpmax*(1/(1 + pow(km/nai,n)))
  ina = 3.0*inapump
  ik = -2.0*inapump
}