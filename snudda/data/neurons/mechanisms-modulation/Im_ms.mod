COMMENT
Mechanism taken from Doron et al., 2017
https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=231427&file=/reproduction/Im.mod#tabs-2

Reference :     Adams et al. 1982 - M-currents and other potassium currents in bullfrog sympathetic neurones

corrected rates using q10 = 2.3, target temperature 34, orginal 21

---------------------------------------------------------------

Neuromodulation is added as functions:
    
    modulationA = 1 + modA*(maxModA-1)*levelA

where:
    
    modA  [0]: is a switch for turning modulation on or off {1/0}
    maxModA [1]: is the maximum modulation for this specific channel (read from the param file)
                    e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}
    levelA  [0]: is an additional parameter for scaling modulation. 
                Can be used simulate non static modulation by gradually changing the value from 0 to 1 {0-1}
									
	  Further neuromodulators can be added by for example:
          modulationA = 1 + modA*(maxModA-1)
	  modulationB = 1 + modB*(maxModB-1)
	  ....

	  etc. for other neuromodulators
	  
	   
								     
[] == default values
{} == ranges
    
ENDCOMMENT

NEURON	{
	SUFFIX Im_ms
	USEION k READ ek WRITE ik
	RANGE gbar, gIm, ik
        RANGE modB, maxModB, levelB
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gbar = 0.00001 (S/cm2) 
	modB = 0
        maxModB = 1 
        levelB = 0
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gIm	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIm = gbar*m*modulationB()
	ik = gIm*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)

	UNITSOFF
		mAlpha = 3.3e-3*exp(2.5*0.04*(v - -35))
		mBeta = 3.3e-3*exp(-2.5*0.04*(v - -35))
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = (1/(mAlpha + mBeta))/qt
	UNITSON
}

FUNCTION modulationB() {
    : returns modulation factor
    
    modulationB = 1 + modB*(maxModB-1)*levelB 
}
