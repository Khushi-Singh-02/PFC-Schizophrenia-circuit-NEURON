NEURON {
  POINT_PROCESS a7
  USEION na WRITE ina
  USEION ca WRITE ica
  USEION k WRITE ik
  POINTER C
  RANGE curr, ica, ina, ik
  RANGE  O1, O2, O3, O4, O5, C5, C4, C3, C2, C1, C0, D0, D1, D2, D3, D4, D5, SD0, SD1, SD2, SD3, SD4, SD5
  RANGE g, gmax, del, tdur
  GLOBAL Erev
  GLOBAL rb1, rb2, rb3, rb4, rb5, Ru1, Ru2, Ru3,Ru4, Ru5, Ro1, Ro2, Ro3, Ro4, Ro5, Rc1, Rc2, Rc3, Rc4, Rc5
  GLOBAL Rd0, Rd1, Rd2, Rd3, Rd4, Rd5, Rr0, Rr1, Rr2, Rr3, Rr4, Rr5, Rd1d0, Rd2d1, Rd3d2, Rd4d3, Rd5d4
  GLOBAL  Rrsd0, Rrsd1, Rrsd2, Rrsd3, Rrsd4, Rrsd5, Rsd0, Rsd1, Rsd2, Rsd3, Rsd4, Rsd5, Rsd1sd0, Rsd2sd1, Rsd3sd2, Rsd4sd3, Rsd5sd4
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
  Erev	= 0   (mV)	: reversal potential
  gmax	= 95  (pS)	: maximal conductance
  tdur = 5 (ms)
  del =  0 (ms)

: Rates

  rb1 = 0.05 (/mM /ms)
  rb2 = 0.04 (/mM /ms)
  rb3 = 0.03 (/mM /ms)
  rb4 = 0.02 (/mM /ms)
  rb5 = 0.01 (/mM /ms)

  Ru1 = 0.5 (/ms)
  Ru2 = 10 (/ms)
  Ru3 = 15 (/ms) 
  Ru4= 20 (/ms)
  Ru5= 25(/ms)

  Rd0 = 0.005 (/ms)
  Rd1 = 0.025 (/ms)
  Rd2 = 0.125 (/ms)
  Rd3 = 0.625 (/ms)
  Rd4 = 3.125 (/ms)
  Rd5 = 15.625 (/ms)

  Rr0 = 5 (/ms)
  Rr1 = 1 (/ms)
  Rr2 = 0.2 (/ms)
  Rr3 = 0.04 (/ms)
  Rr4 = 0.008 (/ms)  
  Rr5 = 0.0016 (/ms)  

  Ro1 = 5e-5 (/ms)
  Ro2 = 1e-3 (/ms)
  Ro3 = 0.02 (/ms)
  Ro4 = 0.4 (/ms)
  Ro5 = 8(/ms)

  Rc1 = 0.02 (/ms)
  Rc2 = 0.02 (/ms)
  Rc3 = 0.02 (/ms)
  Rc4 = 0.02 (/ms)
  Rc5 = 0.02 (/ms)

  Rd1d0 = 0.2 (/ms)
  Rd2d1 = 0.4 (/ms)
  Rd3d2 = 0.6 (/ms)
  Rd4d3 = 0.8 (/ms)
  Rd5d4 = 1 (/ms)

  Rrsd0 = 1e-3 (/ms)
  Rrsd1 = 5e-4 (/ms)
  Rrsd2 = 2.5e-4 (/ms)
  Rrsd3 = 1.25e-4 (/ms)
  Rrsd4 = 6.25e-5 (/ms)
  Rrsd5 = 3.125e-5 (/ms)

  Rsd0 = 2e-5 (/ms)
  Rsd1 = 4e-5 (/ms)
  Rsd2 = 8e-5 (/ms)
  Rsd3 = 16e-5 (/ms)
  Rsd4 = 32e-5 (/ms)
  Rsd5 = 64e-5 (/ms)

  Rsd1sd0 = 0.05 (/ms)
  Rsd2sd1 = 0.1 (/ms)
  Rsd3sd2 = 0.15 (/ms)
  Rsd4sd3 = 0.2 (/ms)
  Rsd5sd4 = 0.25 (/ms)
}

ASSIGNED {
  v		(mV)		: postsynaptic voltage
  curr		(nA)		: current = g*(v - Erev)
  g 		(pS)		: conductance
  tr (mM)				:ACh conc
  C (mM)				:pointer
  Rb1		(/ms)    
  Rb2		(/ms)    
  Rb3		(/ms)    
  Rb4		(/ms)    
  Rb5		(/ms)    
  ica (nA)
  ina (nA)
  ik (nA)
}

STATE {
  C5 C4 C3 C2 C1 C0 
  O1 O2 O3 O4 O5 
  D0 D1 D2 D3 D4 D5 
  SD0 SD1 SD2 SD3 SD4 SD5
}

INITIAL {
  C5 =1
  C4 =0
  C3 =0 
  C2 =0 
  C1 =0 
  C0 =0 
  O1 =0 
  O2 =0 
  O3 =0 
  O4 =0 
  O5 =0 
  D0 =0 
  D1 =0 
  D2 =0 
  D3 =0 
  D4 =0 
  D5 =0 
  SD0 =0 
  SD1 =0 
  SD2 =0 
  SD3 =0 
  SD4 =0 
  SD5 =0 
  tr = 0
}

BREAKPOINT {
  SOLVE kstates METHOD sparse

  g = gmax * (O1 + O2 + O3 + O4 + O5)
  curr = (1e-6) * g * (v - Erev)  
  ica = 0.77*curr 
  ina = 0.115*curr
  ik = 0.115*curr
}

KINETIC kstates {
    Rb1 = rb1 * C
    Rb2 = rb2 * C
    Rb3 = rb3 * C
    Rb4 = rb4 * C
    Rb5 = rb5 * C
    
   ~ C5 <-> C4 (Rb1, Ru1 )
   ~ C4 <-> C3 (Rb2, Ru2 )
   ~ C3 <-> C2 (Rb3, Ru3 )
   ~ C2 <-> C1 (Rb4, Ru4 )
   ~ C1 <-> C0 (Rb5, Ru5 )

   ~ C4 <-> O1 (Ro1, Rc1 )
   ~ C3 <-> O2 (Ro2, Rc2 )
   ~ C2 <-> O3 (Ro3, Rc3 )
   ~ C1 <-> O4 (Ro4, Rc4 )
   ~ C0 <-> O5 (Ro5, Rc5 )

   ~ C5 <-> D0 (Rd0, Rr0 )
   ~ C4 <-> D1 (Rd1, Rr1 )
   ~ C3 <-> D2 (Rd2, Rr2 )
   ~ C2 <-> D3 (Rd3, Rr3 )
   ~ C1 <-> D4 (Rd4, Rr4 )
   ~ C0 <-> D5 (Rd5, Rr5 )

   ~ D0 <-> D1 (Rb1, Rd1d0)
   ~ D1 <-> D2 (Rb2, Rd2d1)
   ~ D2 <-> D3 (Rb1, Rd3d2)
   ~ D0 <-> D1 (Rb1, Rd4d3)
   ~ D0 <-> D1 (Rb1, Rd5d4)

   ~ D0 <-> SD0 (Rsd0, Rrsd0)
   ~ D1 <-> SD1 (Rsd1, Rrsd1)
   ~ D2 <-> SD2 (Rsd2, Rrsd2)
   ~ D3 <-> SD3 (Rsd3, Rrsd3)
   ~ D4 <-> SD4 (Rsd4, Rrsd4)
   ~ D5 <-> SD5 (Rsd5, Rrsd5)

   ~ SD0 <-> SD1 (Rb1, Rsd1sd0)
   ~ SD1 <-> SD2 (Rb2, Rsd2sd1)
   ~ SD2 <-> SD3 (Rb3, Rsd3sd2)
   ~ SD3 <-> SD4 (Rb4, Rsd4sd3)
   ~ SD4 <-> SD5 (Rb5, Rsd5sd4)

  CONSERVE C5 + C4 + C3 + C2 + C1 + C0 + O1 + O2 + O3 + O4 + O5 + D0 + D1 + D2 + D3 + D4 + D5 + SD0 + SD1 + SD2 + SD3 + SD4 +SD5= 1
 }