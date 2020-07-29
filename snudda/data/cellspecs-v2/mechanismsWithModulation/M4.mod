TITLE Mod file for component: Component(id=M4)

COMMENT 
	Subcellular cascade automatically generated from SBML.
		Generated using sbml2mod_converter.py and manually corrected for concentrations
		date: 2019-09-20 
ENDCOMMENT

NEURON { 
	THREADSAFE
	POINT_PROCESS M4 
	RANGE Ca,tau1, tau2, Ach
        POINTER conc_ACH

}

PARAMETER {

	
	kf_AchxM4R            =  3.38e-4 
	kb_AchxM4R            =  2.5
	kf_AchxM4R_Gi         =  0.036
	kbAchxM4R_Gi          =  20.0
	kf_M4RxGiabg          =  0.00178
	kb_M4RxGiabg          =  0.42
	kf_AchM4RxGiabg       =  0.192
	kb_AchM4RxGiabg       =  3.36
	kf_Ach_M4R_Gi         =  8.0
	kb_Ach_M4R_Gi         =  0
	kf_GiaGTP             =  4.0
	kb_GiaGTP             =  0
	kf_GiaGDP             =  2.0
	kb_GiaGDP             =  0
	init_AChE	      =  1 
	init_Ach              =  105 
	init_M4R              =  120 
	init_Ach_M4R          =  0 
	init_Giabg            =  1680 
	init_Giabg_M4R        =  80 
	init_Ach_M4R_Giabg    =  0 
	init_GiaGTP           =  80 
	init_GiaGDP           =  0
	init_Gbg              =  120 

	: --------- concentration converter (cc): nM -> mM; time converter s -> ms
	cc       = 1e-6      
	two_s_rs = 1e3  (kHz)
	one_s_rs = 1e-3 (kHz)
	tau1 = 100 (ms) <1e-9,1e9>
	tau2 = 500 (ms) <1e-9,1e9>
	e=0	(mV)

}

STATE {

	Ach
	M4R
	Ach_M4R
	Giabg
	Giabg_M4R
	Ach_M4R_Giabg
	GiaGTP
	GiaGDP
	Gbg
	A
	B

}
ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	conc_ACH (1)
}


BREAKPOINT { 	
	SOLVE kin METHOD sparse
	total_M4R()
	Total_Gi()
	
	 
	Ach = acetylcholine()
				
}

INITIAL {

	LOCAL tp
	if (tau1/tau2 > 0.9999) {
		tau1 = 0.9999*tau2
	}
	if (tau1/tau2 < 1e-9) {
		tau1 = tau2*1e-9
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor

	Ach                   = init_Ach * cc
	M4R                   = init_M4R * cc
	Ach_M4R               = init_Ach_M4R * cc
	Giabg                 = init_Giabg * cc
	Giabg_M4R             = init_Giabg_M4R * cc
	Ach_M4R_Giabg         = init_Ach_M4R_Giabg * cc
	GiaGTP                = init_GiaGTP * cc
	GiaGDP                = init_GiaGDP * cc
	Gbg                   = init_Gbg * cc

	

}

KINETIC kin {

	~ Ach + M4R           <->  Ach_M4R                   		(kf_AchxM4R*two_s_rs,         kb_AchxM4R*one_s_rs)
	~ Ach_M4R + Giabg     <->  Ach_M4R_Giabg             		(kf_AchM4RxGiabg*two_s_rs,    kb_AchM4RxGiabg*one_s_rs)
	~ Ach_M4R_Giabg       <->  GiaGTP + Ach_M4R + Gbg    		(kf_Ach_M4R_Gi*one_s_rs,      kb_Ach_M4R_Gi)
	~ M4R + Giabg         <->  Giabg_M4R                 		(kf_M4RxGiabg*two_s_rs,       kb_M4RxGiabg*one_s_rs)
	~ GiaGTP              <->  GiaGDP                    		(kf_GiaGTP*one_s_rs,          kb_GiaGTP)
	~ GiaGDP + Gbg        <->  Giabg                     		(kf_GiaGDP*one_s_rs,          kb_GiaGDP)
	~ Ach + Giabg_M4R     <->  Ach_M4R_Giabg             		(kf_AchxM4R_Gi*two_s_rs,      kbAchxM4R_Gi*one_s_rs)

}

FUNCTION total_M4R() {
	total_M4R = M4R + Ach_M4R + Ach_M4R_Giabg + Giabg_M4R
}

FUNCTION Total_Gi() {
	Total_Gi = Giabg + Giabg_M4R + GiaGTP + GiaGDP + Ach_M4R_Giabg
}

FUNCTION acetylcholine() {
     acetylcholine = conc_ACH
}
