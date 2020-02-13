TITLE Calcium dynamics for L and T calcium pool

NEURON {
    SUFFIX caldyn_ms
    USEION cal READ ical, cali WRITE cali VALENCE 2
    RANGE pump, cainf, taur, drive, depth
}

UNITS {
    (molar) = (1/liter) 
    (mM) = (millimolar)
    (um) = (micron)
    (mA) = (milliamp)
    (msM) = (ms mM)
    FARADAY = (faraday) (coulomb)
}

PARAMETER {
    drive = 10000 (1)
    depth = 0.2  (um)
    cainf = 70e-6 (mM)
    taur = 43 (ms)
    kt = 1e-4 (mM/ms)
    kd = 1e-4 (mM)
    pump = 0.02
}

STATE { cali (mM) }

INITIAL { cali = cainf }

ASSIGNED {
    ical (mA/cm2)
    drive_channel (mM/ms)
    drive_pump (mM/ms)
}
    
BREAKPOINT {
    SOLVE state METHOD cnexp
}


DERIVATIVE state { 
    
    : force concentration to stay above cainf by only pumping if larger
    drive_channel = -drive*ical/(2*FARADAY*depth)
    drive_pump    = -kt*(cali-cainf)/(cali+kd)
    
    if (drive_channel <= 0.) { drive_channel = 0. }
    
    cali' = drive_channel + pump*drive_pump + (cainf-cali)/taur
}

COMMENT

Original NEURON model by Wolf (2005) and Destexhe (1992).  Adaptation by
Alexander Kozlov <akozlov@kth.se>. Updated by Robert Lindroos <robert.lindroos@ki.se>.

Updates by RL:
-cainf changed from 10 to 70 nM (sabatini et al., 2002 The Life Cycle of Ca 2+ Ions in Dendritic Spines)
-pump updated to only be active if cai > cainf (neutralized by adding reversed entity) 

[1] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,
O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions
and entrainment to oscillations in a computational model of the nucleus
accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.

ENDCOMMENT
