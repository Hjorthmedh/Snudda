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
 
#define nrn_init _nrn_init__Kv3_ch
#define _nrn_initial _nrn_initial__Kv3_ch
#define nrn_cur _nrn_cur__Kv3_ch
#define _nrn_current _nrn_current__Kv3_ch
#define nrn_jacob _nrn_jacob__Kv3_ch
#define nrn_state _nrn_state__Kv3_ch
#define _net_receive _net_receive__Kv3_ch 
#define _f_settables _f_settables__Kv3_ch 
#define settables settables__Kv3_ch 
#define states states__Kv3_ch 
 
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
#define gmax _p[0]
#define iKv3 _p[1]
#define m _p[2]
#define h _p[3]
#define ek _p[4]
#define Dm _p[5]
#define Dh _p[6]
#define ik _p[7]
#define minf _p[8]
#define taum _p[9]
#define hinf _p[10]
#define tauh _p[11]
#define v _p[12]
#define _g _p[13]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
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
 static void _hoc_settables(void);
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
 "setdata_Kv3_ch", _hoc_setdata,
 "settables_Kv3_ch", _hoc_settables,
 0, 0
};
 
static void _check_settables(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_settables(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
#define h0 h0_Kv3_ch
 double h0 = 0.6;
#define k_h k_h_Kv3_ch
 double k_h = -10;
#define k_m k_m_Kv3_ch
 double k_m = 7.8;
#define phi_h phi_h_Kv3_ch
 double phi_h = 0;
#define phi_m phi_m_Kv3_ch
 double phi_m = -26;
#define sigma_h1 sigma_h1_Kv3_ch
 double sigma_h1 = -10;
#define sigma_h0 sigma_h0_Kv3_ch
 double sigma_h0 = 10;
#define sigma_m1 sigma_m1_Kv3_ch
 double sigma_m1 = -12;
#define sigma_m0 sigma_m0_Kv3_ch
 double sigma_m0 = 13;
#define tau_h1 tau_h1_Kv3_ch
 double tau_h1 = 33;
#define tau_h0 tau_h0_Kv3_ch
 double tau_h0 = 7;
#define theta_h theta_h_Kv3_ch
 double theta_h = -20;
#define tau_m1 tau_m1_Kv3_ch
 double tau_m1 = 14;
#define tau_m0 tau_m0_Kv3_ch
 double tau_m0 = 0.1;
#define theta_m theta_m_Kv3_ch
 double theta_m = -26;
#define usetable usetable_Kv3_ch
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_Kv3_ch", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "theta_m_Kv3_ch", "mV",
 "k_m_Kv3_ch", "mV",
 "tau_m0_Kv3_ch", "ms",
 "tau_m1_Kv3_ch", "ms",
 "phi_m_Kv3_ch", "mV",
 "sigma_m0_Kv3_ch", "mV",
 "sigma_m1_Kv3_ch", "mV",
 "theta_h_Kv3_ch", "mV",
 "k_h_Kv3_ch", "mV",
 "tau_h0_Kv3_ch", "ms",
 "tau_h1_Kv3_ch", "ms",
 "phi_h_Kv3_ch", "mV",
 "sigma_h0_Kv3_ch", "mV",
 "sigma_h1_Kv3_ch", "mV",
 "gmax_Kv3_ch", "mho/cm2",
 "iKv3_Kv3_ch", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double m0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "theta_m_Kv3_ch", &theta_m_Kv3_ch,
 "k_m_Kv3_ch", &k_m_Kv3_ch,
 "tau_m0_Kv3_ch", &tau_m0_Kv3_ch,
 "tau_m1_Kv3_ch", &tau_m1_Kv3_ch,
 "phi_m_Kv3_ch", &phi_m_Kv3_ch,
 "sigma_m0_Kv3_ch", &sigma_m0_Kv3_ch,
 "sigma_m1_Kv3_ch", &sigma_m1_Kv3_ch,
 "h0_Kv3_ch", &h0_Kv3_ch,
 "theta_h_Kv3_ch", &theta_h_Kv3_ch,
 "k_h_Kv3_ch", &k_h_Kv3_ch,
 "tau_h0_Kv3_ch", &tau_h0_Kv3_ch,
 "tau_h1_Kv3_ch", &tau_h1_Kv3_ch,
 "phi_h_Kv3_ch", &phi_h_Kv3_ch,
 "sigma_h0_Kv3_ch", &sigma_h0_Kv3_ch,
 "sigma_h1_Kv3_ch", &sigma_h1_Kv3_ch,
 "usetable_Kv3_ch", &usetable_Kv3_ch,
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
"Kv3_ch",
 "gmax_Kv3_ch",
 "iKv3_Kv3_ch",
 0,
 0,
 "m_Kv3_ch",
 "h_Kv3_ch",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 14, _prop);
 	/*initialize range parameters*/
 	gmax = 0.001;
 	iKv3 = 0;
 	_prop->param = _p;
 	_prop->param_size = 14;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
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

 void _Kv3_ch_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("k", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 14, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 Kv3_ch /home/akozlov/doc/hbp-bsp-live-papers-dev-priv/2020/hjorth_et_al_2020/work/Snudda/x86_64/Kv3_ch.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_minf;
 static double *_t_taum;
 static double *_t_hinf;
 static double *_t_tauh;
static int _reset;
static char *modelname = "fast activated potassium Kv3 (Kv3.1/3.4) channel for GPe neuron";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_settables(_threadargsprotocomma_ double);
static int settables(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_settables(_threadargsprotocomma_ double _lv);
 static int _slist1[2], _dlist1[2];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   settables ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / taum ;
   Dh = ( hinf - h ) / tauh ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 settables ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taum )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tauh )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   settables ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taum)))*(- ( ( ( minf ) ) / taum ) / ( ( ( ( - 1.0 ) ) ) / taum ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tauh)))*(- ( ( ( hinf ) ) / tauh ) / ( ( ( ( - 1.0 ) ) ) / tauh ) - h) ;
   }
  return 0;
}
 static double _mfac_settables, _tmin_settables;
  static void _check_settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_settables =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_settables)/400.; _mfac_settables = 1./_dx;
   for (_i=0, _x=_tmin_settables; _i < 401; _x += _dx, _i++) {
    _f_settables(_p, _ppvar, _thread, _nt, _x);
    _t_minf[_i] = minf;
    _t_taum[_i] = taum;
    _t_hinf[_i] = hinf;
    _t_tauh[_i] = tauh;
   }
  }
 }

 static int settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_settables(_p, _ppvar, _thread, _nt);
#endif
 _n_settables(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_settables(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_settables(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_settables * (_lv - _tmin_settables);
 if (isnan(_xi)) {
  minf = _xi;
  taum = _xi;
  hinf = _xi;
  tauh = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 taum = _t_taum[0];
 hinf = _t_hinf[0];
 tauh = _t_tauh[0];
 return; }
 if (_xi >= 400.) {
 minf = _t_minf[400];
 taum = _t_taum[400];
 hinf = _t_hinf[400];
 tauh = _t_tauh[400];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 taum = _t_taum[_i] + _theta*(_t_taum[_i+1] - _t_taum[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 tauh = _t_tauh[_i] + _theta*(_t_tauh[_i+1] - _t_tauh[_i]);
 }

 
static int  _f_settables ( _threadargsprotocomma_ double _lv ) {
   minf = 1.0 / ( 1.0 + exp ( ( theta_m - _lv ) / k_m ) ) ;
   taum = tau_m0 + ( tau_m1 - tau_m0 ) / ( exp ( ( phi_m - _lv ) / sigma_m0 ) + exp ( ( phi_m - _lv ) / sigma_m1 ) ) ;
   hinf = h0 + ( 1.0 - h0 ) / ( 1.0 + exp ( ( theta_h - _lv ) / k_h ) ) ;
   tauh = tau_h0 + ( tau_h1 - tau_h0 ) / ( exp ( ( phi_h - _lv ) / sigma_h0 ) + exp ( ( phi_h - _lv ) / sigma_h1 ) ) ;
    return 0; }
 
static void _hoc_settables(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 settables ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 2;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 2; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
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
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
 {
   settables ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
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

#if 0
 _check_settables(_p, _ppvar, _thread, _nt);
#endif
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
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ik = gmax * m * m * m * m * h * ( v - ek ) ;
   iKv3 = ik ;
   }
 _current += ik;

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
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
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
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  } }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
   _t_minf = makevector(401*sizeof(double));
   _t_taum = makevector(401*sizeof(double));
   _t_hinf = makevector(401*sizeof(double));
   _t_tauh = makevector(401*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/akozlov/doc/hbp-bsp-live-papers-dev-priv/2020/hjorth_et_al_2020/work/Snudda/snudda/data/cellspecs/mechanisms/Kv3_ch.mod";
static const char* nmodl_file_text = 
  "TITLE fast activated potassium Kv3 (Kv3.1/3.4) channel for GPe neuron\n"
  "\n"
  "COMMENT\n"
  " modeled by Gunay et al., 2008\n"
  " implemented in NEURON by Kitano, 2011\n"
  "ENDCOMMENT\n"
  "\n"
  "UNITS {\n"
  " (mV) = (millivolt)\n"
  " (mA) = (milliamp)\n"
  "}\n"
  "\n"
  "NEURON {\n"
  " SUFFIX Kv3_ch\n"
  " USEION k READ ek WRITE ik\n"
  " RANGE gmax, iKv3\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  " v (mV)\n"
  " dt (ms)\n"
  " gmax  = 0.001 (mho/cm2)\n"
  " iKv3  = 0.0 (mA/cm2)\n"
  " ek (mV)\n"
  "\n"
  " theta_m = -26.0 (mV)\n"
  " k_m = 7.8 (mV)\n"
  " tau_m0 = 0.1 (ms)\n"
  " tau_m1 = 14.0 (ms)\n"
  " phi_m = -26.0 (mV)\n"
  " sigma_m0 = 13.0 (mV)\n"
  " sigma_m1 = -12.0 (mV)\n"
  "\n"
  " h0 = 0.6\n"
  " theta_h = -20.0 (mV)\n"
  " k_h = -10.0 (mV)\n"
  " tau_h0 = 7.0 (ms)\n"
  " tau_h1 = 33.0 (ms)\n"
  " phi_h = 0.0 (mV)\n"
  " sigma_h0 = 10.0 (mV)\n"
  " sigma_h1 = -10.0 (mV)\n"
  "}\n"
  "\n"
  "STATE {\n"
  " m h\n"
  "}\n"
  "\n"
  "ASSIGNED { \n"
  " ik (mA/cm2)\n"
  " minf\n"
  " taum (ms)\n"
  " hinf\n"
  " tauh (ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  " SOLVE states METHOD cnexp\n"
  " ik  = gmax*m*m*m*m*h*(v-ek)\n"
  " iKv3 = ik\n"
  "}\n"
  "\n"
  "UNITSOFF\n"
  "\n"
  "INITIAL {\n"
  " settables(v)\n"
  " m = minf\n"
  " h = hinf\n"
  "}\n"
  "\n"
  "DERIVATIVE states {  \n"
  " settables(v)\n"
  " m' = (minf - m)/taum\n"
  " h' = (hinf - h)/tauh\n"
  "}\n"
  "\n"
  "PROCEDURE settables(v) {\n"
  "        TABLE minf, taum, hinf, tauh FROM -100 TO 100 WITH 400\n"
  "\n"
  "	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))\n"
  "	taum = tau_m0 + (tau_m1 - tau_m0)/(exp((phi_m - v)/sigma_m0) + exp((phi_m - v)/sigma_m1))\n"
  "	hinf = h0 + (1.0 - h0) / (1.0 + exp((theta_h - v)/k_h))\n"
  "	tauh = tau_h0 + (tau_h1 - tau_h0)/(exp((phi_h - v)/sigma_h0) + exp((phi_h - v)/sigma_h1))\n"
  "}\n"
  "\n"
  "UNITSON\n"
  ;
#endif
