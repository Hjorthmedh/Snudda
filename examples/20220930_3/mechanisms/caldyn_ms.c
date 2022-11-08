/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__caldyn_ms
#define _nrn_initial _nrn_initial__caldyn_ms
#define nrn_cur _nrn_cur__caldyn_ms
#define _nrn_current _nrn_current__caldyn_ms
#define nrn_jacob _nrn_jacob__caldyn_ms
#define nrn_state _nrn_state__caldyn_ms
#define _net_receive _net_receive__caldyn_ms 
#define state state__caldyn_ms 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define drive _p[0]
#define depth _p[1]
#define cainf _p[2]
#define taur _p[3]
#define pump _p[4]
#define cali _p[5]
#define Dcali _p[6]
#define ical _p[7]
#define drive_channel _p[8]
#define drive_pump _p[9]
#define v _p[10]
#define _g _p[11]
#define _ion_ical	*_ppvar[0]._pval
#define _ion_cali	*_ppvar[1]._pval
#define _style_cal	*((int*)_ppvar[2]._pvoid)
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_caldyn_ms", _hoc_setdata,
 0, 0
};
 /* declare global and static user variables */
#define kd kd_caldyn_ms
 double kd = 0.0001;
#define kt kt_caldyn_ms
 double kt = 0.0001;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "kt_caldyn_ms", "mM/ms",
 "kd_caldyn_ms", "mM",
 "drive_caldyn_ms", "1",
 "depth_caldyn_ms", "um",
 "cainf_caldyn_ms", "mM",
 "taur_caldyn_ms", "ms",
 0,0
};
 static double cali0 = 0;
 static double delta_t = 0.01;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "kt_caldyn_ms", &kt_caldyn_ms,
 "kd_caldyn_ms", &kd_caldyn_ms,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[3]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"caldyn_ms",
 "drive_caldyn_ms",
 "depth_caldyn_ms",
 "cainf_caldyn_ms",
 "taur_caldyn_ms",
 "pump_caldyn_ms",
 0,
 0,
 0,
 0};
 static Symbol* _cal_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	drive = 10000;
 	depth = 0.2;
 	cainf = 7e-05;
 	taur = 43;
 	pump = 0.02;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_cal_sym);
 nrn_check_conc_write(_prop, prop_ion, 1);
 nrn_promote(prop_ion, 3, 0);
 	_ppvar[0]._pval = &prop_ion->param[3]; /* ical */
 	_ppvar[1]._pval = &prop_ion->param[1]; /* cali */
 	_ppvar[2]._pvoid = (void*)(&(prop_ion->dparam[0]._i)); /* iontype for cal */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _caldyn_ms_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("cal", 2.0);
 	_cal_sym = hoc_lookup("cal_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "cal_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "#cal_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	nrn_writes_conc(_mechtype, 0);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 caldyn_ms C:/Users/bo.bekkouche/PycharmProjects/currentinjection/Snudda/examples/bgd01/parkinson/20211105/PD0/mechanisms/caldyn_ms.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0xb.c72aa8304416p+13, 96485.3}; /* 96485.3321233100141 */
static int _reset;
static char *modelname = "Calcium dynamics for L and T calcium pool";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   drive_channel = - drive * ical / ( 2.0 * FARADAY * depth ) ;
   drive_pump = - kt * ( cali - cainf ) / ( cali + kd ) ;
   if ( drive_channel <= 0. ) {
     drive_channel = 0. ;
     }
   Dcali = drive_channel + pump * drive_pump + ( cainf - cali ) / taur ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 drive_channel = - drive * ical / ( 2.0 * FARADAY * depth ) ;
 drive_pump = - kt * ( cali - cainf ) / ( cali + kd ) ;
 if ( drive_channel <= 0. ) {
   drive_channel = 0. ;
   }
 Dcali = Dcali  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taur )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   drive_channel = - drive * ical / ( 2.0 * FARADAY * depth ) ;
   drive_pump = - kt * ( cali - cainf ) / ( cali + kd ) ;
   if ( drive_channel <= 0. ) {
     drive_channel = 0. ;
     }
    cali = cali + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taur)))*(- ( drive_channel + ( pump )*( drive_pump ) + ( ( cainf ) ) / taur ) / ( ( ( ( - 1.0 ) ) ) / taur ) - cali) ;
   }
  return 0;
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ical = _ion_ical;
  cali = _ion_cali;
  cali = _ion_cali;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  _ion_cali = cali;
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 	_pv[0] = &(_ion_cali);
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ical = _ion_ical;
  cali = _ion_cali;
  cali = _ion_cali;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_cal_sym, _ppvar, 0, 3);
   nrn_update_ion_pointer(_cal_sym, _ppvar, 1, 1);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
 {
   cali = cainf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ical = _ion_ical;
  cali = _ion_cali;
  cali = _ion_cali;
 initmodel(_p, _ppvar, _thread, _nt);
  _ion_cali = cali;
  nrn_wrote_conc(_cal_sym, (&(_ion_cali)) - 1, _style_cal);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{
} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 
}
 
}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ical = _ion_ical;
  cali = _ion_cali;
  cali = _ion_cali;
 {   state(_p, _ppvar, _thread, _nt);
  } {
   }
  _ion_cali = cali;
}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(cali) - _p;  _dlist1[0] = &(Dcali) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "caldyn_ms.mod";
static const char* nmodl_file_text = 
  "TITLE Calcium dynamics for L and T calcium pool\n"
  "\n"
  "NEURON {\n"
  "    SUFFIX caldyn_ms\n"
  "    USEION cal READ ical, cali WRITE cali VALENCE 2\n"
  "    RANGE pump, cainf, taur, drive, depth\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "    (molar) = (1/liter) \n"
  "    (mM) = (millimolar)\n"
  "    (um) = (micron)\n"
  "    (mA) = (milliamp)\n"
  "    (msM) = (ms mM)\n"
  "    FARADAY = (faraday) (coulomb)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "    drive = 10000 (1)\n"
  "    depth = 0.2  (um)\n"
  "    cainf = 70e-6 (mM)\n"
  "    taur = 43 (ms)\n"
  "    kt = 1e-4 (mM/ms)\n"
  "    kd = 1e-4 (mM)\n"
  "    pump = 0.02\n"
  "}\n"
  "\n"
  "STATE { cali (mM) }\n"
  "\n"
  "INITIAL { cali = cainf }\n"
  "\n"
  "ASSIGNED {\n"
  "    ical (mA/cm2)\n"
  "    drive_channel (mM/ms)\n"
  "    drive_pump (mM/ms)\n"
  "}\n"
  "    \n"
  "BREAKPOINT {\n"
  "    SOLVE state METHOD cnexp\n"
  "}\n"
  "\n"
  "\n"
  "DERIVATIVE state { \n"
  "    \n"
  "    : force concentration to stay above cainf by only pumping if larger\n"
  "    drive_channel = -drive*ical/(2*FARADAY*depth)\n"
  "    drive_pump    = -kt*(cali-cainf)/(cali+kd)\n"
  "    \n"
  "    if (drive_channel <= 0.) { drive_channel = 0. }\n"
  "    \n"
  "    cali' = drive_channel + pump*drive_pump + (cainf-cali)/taur\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Original NEURON model by Wolf (2005) and Destexhe (1992).  Adaptation by\n"
  "Alexander Kozlov <akozlov@kth.se>. Updated by Robert Lindroos <robert.lindroos@ki.se>.\n"
  "\n"
  "Updates by RL:\n"
  "-cainf changed from 10 to 70 nM (sabatini et al., 2002 The Life Cycle of Ca 2+ Ions in Dendritic Spines)\n"
  "-pump updated to only be active if cai > cainf (neutralized by adding reversed entity) \n"
  "\n"
  "[1] Wolf JA, Moyer JT, Lazarewicz MT, Contreras D, Benoit-Marand M,\n"
  "O'Donnell P, Finkel LH (2005) NMDA/AMPA ratio impacts state transitions\n"
  "and entrainment to oscillations in a computational model of the nucleus\n"
  "accumbens medium spiny projection neuron. J Neurosci 25(40):9080-95.\n"
  "\n"
  "ENDCOMMENT\n"
  ;
#endif
