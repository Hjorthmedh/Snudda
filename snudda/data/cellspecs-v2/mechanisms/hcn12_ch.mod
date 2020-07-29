COMMENT

Josh Held's adaptation to suit HCN1+2.  12/22/2003

****
Kinetic model of HCN2 channel gating from Wang et al 2002.

In this model channel opening is coupled to a change in the affinity of the cyclic nucleotide binding domain for cAMP which is manifest as a shift in the activation curve toward more positive potentials.  This model explains the slow activation kinetics of Ih associated with low concentrations of cAMP.

For further details email Matt Nolan at mfnolan@fido.cpmc.columbia.edu.

Reference

Wang J., Chen S., Nolan M.F. and Siegelbaum S.A. (2002). Activity-dependent regulation of HCN pacemaker channels by cyclicAMP: signalling through dynamic allosteric coupling. Neuron 36, 1-20.
****

neuromodulation is added as functions:
    
    modulation = 1 + damod*(maxMod-1)

where:
    
    damod  [0]: is a switch for turning modulation on or off {1/0}
    maxMod [1]: is the maximum modulation for this specific channel (read from the param file)
                    e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}

[] == default values
{} == ranges
    
ENDCOMMENT

NEURON {
	SUFFIX hcn12_ch
	NONSPECIFIC_CURRENT i
	RANGE i, ehcn, g, gbar
	GLOBAL a0, b0, ah, bh, ac, bc, aa0, ba0
	GLOBAL aa0, ba0, aah, bah, aac, bac
	GLOBAL kon, koff, b, bf, ai, gca, shift
    RANGE damod, maxMod
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	gbar    = 1		(S/cm2)
	ehcn    = -20		(mV)
	a0      = .006		(/ms)	: parameters for alpha and beta
	b0      = .0008		(/ms)
	ah      = -96		(mV)
	bh      = -51.7		(mV)
	ac      = -.155		(/mV)
	bc      = .144		(/mV)
	aa0     = .0006		(/ms)	: parameters for alphaa and betaa
	ba0     = .004		(/ms)
	aah     = -94.2		(mV)
	bah     = -35.5		(mV)
	aac     = -.075		(/mV)
	bac     = .144		(/mV)
	kon     = 30		(/mM-ms) : cyclic AMP binding parameters
	koff    = 4.5e-05	(/ms)
	b       = 80
	bf      = 8.94
	ai	= 1e-05		(mM)	:concentration cyclic AMP
	gca     = 1			: relative conductance of the bound state
	shift   = -17		(mV)	: shift in voltage dependence
	q10v    = 4                     : q10 value from Magee 1998
	q10a    = 1.5			: estimated q10 for the cAMP binding reaction
	celsius			(degC)
    damod = 0
    maxMod = 1
}

ASSIGNED {
	v	(mV)
	g	(S/cm2)
	i	(mA/cm2)
	alpha	(/ms)
	beta    (/ms)
	alphaa	(/ms)
	betaa	(/ms)
}

STATE {
	c
	cac
	o
	cao
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	g = gbar*(o + cao*gca)*modulation()
	i = g*(v-ehcn)
}

KINETIC kin {
	LOCAL qa
	qa = q10a^((celsius-22 (degC))/10 (degC))
	rates(v)
	~ c <-> o       (alpha, beta)
	~ c <-> cac     (kon*qa*ai/bf,koff*qa*b/bf)
	~ o <-> cao     (kon*qa*ai, koff*qa)
	~ cac <-> cao   (alphaa, betaa)
	CONSERVE c + cac + o + cao = 1
}

PROCEDURE rates(v(mV)) {
	LOCAL qv
	qv = q10v^((celsius-22 (degC))/10 (degC))
	if (v > -200) {
		alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))
	} else {
		alpha = a0*qv / (1 + exp(-((-200)-ah-shift)*ac))
		beta = b0*qv / (1 + exp(-((-200)-bh-shift)*bc))
		alphaa = aa0*qv / (1 + exp(-((-200)-aah-shift)*aac))
		betaa = ba0*qv / (1 + exp(-((-200)-bah-shift)*bac))
	}
}

FUNCTION modulation() {
    : returns modulation factor
    
    modulation = 1 + damod*(maxMod-1)
}
