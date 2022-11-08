/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
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
 
#define nrn_init _nrn_init__bk_ch
#define _nrn_initial _nrn_initial__bk_ch
#define nrn_cur _nrn_cur__bk_ch
#define _nrn_current _nrn_current__bk_ch
#define nrn_jacob _nrn_jacob__bk_ch
#define nrn_state _nrn_state__bk_ch
#define _net_receive _net_receive__bk_ch 
#define rates rates__bk_ch 
#define states states__bk_ch 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gbar _p[0]
#define g _p[1]
#define ik _p[2]
#define m _p[3]
#define z _p[4]
#define h _p[5]
#define ek _p[6]
#define cai _p[7]
#define Dm _p[8]
#define Dz _p[9]
#define Dh _p[10]
#define _g _p[11]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
#define _ion_cai	*_ppvar[3]._pval
 
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
 /* external NEURON variables */
 /* declaration of user functions */
 static void _hoc_rates(void);
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
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_bk_ch", _hoc_setdata,
 "rates_bk_ch", _hoc_rates,
 0, 0
};
 /* declare global and static user variables */
#define htau_k2 htau_k2_bk_ch
 double htau_k2 = 5.2;
#define htau_k1 htau_k1_bk_ch
 double htau_k1 = -12.9;
#define htau_vh2 htau_vh2_bk_ch
 double htau_vh2 = 48.5;
#define htau_vh1 htau_vh1_bk_ch
 double htau_vh1 = -54.2;
#define htau_y0 htau_y0_bk_ch
 double htau_y0 = 0.0019;
#define h_k h_k_bk_ch
 double h_k = 5.8;
#define h_vh h_vh_bk_ch
 double h_vh = -32;
#define h_y0 h_y0_bk_ch
 double h_y0 = 0.085;
#define htau htau_bk_ch
 double htau = 0;
#define hinf hinf_bk_ch
 double hinf = 0;
#define mtau_k2 mtau_k2_bk_ch
 double mtau_k2 = 10.1;
#define mtau_k1 mtau_k1_bk_ch
 double mtau_k1 = -10;
#define mtau_vh2 mtau_vh2_bk_ch
 double mtau_vh2 = 86.4;
#define mtau_vh1 mtau_vh1_bk_ch
 double mtau_vh1 = -33.3;
#define mtau_y0 mtau_y0_bk_ch
 double mtau_y0 = 0.000505;
#define m_k m_k_bk_ch
 double m_k = 6.2;
#define m_vh m_vh_bk_ch
 double m_vh = -28.9;
#define mtau mtau_bk_ch
 double mtau = 0;
#define minf minf_bk_ch
 double minf = 0;
#define z_coef z_coef_bk_ch
 double z_coef = 0.001;
#define ztau ztau_bk_ch
 double ztau = 1;
#define zinf zinf_bk_ch
 double zinf = 0;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "m_vh_bk_ch", "mV",
 "m_k_bk_ch", "mV",
 "mtau_y0_bk_ch", "s",
 "mtau_vh1_bk_ch", "mV",
 "mtau_k1_bk_ch", "mV",
 "mtau_vh2_bk_ch", "mV",
 "mtau_k2_bk_ch", "mV",
 "z_coef_bk_ch", "mM",
 "ztau_bk_ch", "ms",
 "h_vh_bk_ch", "mV",
 "h_k_bk_ch", "mV",
 "htau_y0_bk_ch", "s",
 "htau_vh1_bk_ch", "mV",
 "htau_k1_bk_ch", "mV",
 "htau_vh2_bk_ch", "mV",
 "htau_k2_bk_ch", "mV",
 "mtau_bk_ch", "ms",
 "htau_bk_ch", "ms",
 "gbar_bk_ch", "mho/cm2",
 "g_bk_ch", "S/cm2",
 "ik_bk_ch", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double v = 0;
 static double z0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "m_vh_bk_ch", &m_vh_bk_ch,
 "m_k_bk_ch", &m_k_bk_ch,
 "mtau_y0_bk_ch", &mtau_y0_bk_ch,
 "mtau_vh1_bk_ch", &mtau_vh1_bk_ch,
 "mtau_k1_bk_ch", &mtau_k1_bk_ch,
 "mtau_vh2_bk_ch", &mtau_vh2_bk_ch,
 "mtau_k2_bk_ch", &mtau_k2_bk_ch,
 "z_coef_bk_ch", &z_coef_bk_ch,
 "ztau_bk_ch", &ztau_bk_ch,
 "h_y0_bk_ch", &h_y0_bk_ch,
 "h_vh_bk_ch", &h_vh_bk_ch,
 "h_k_bk_ch", &h_k_bk_ch,
 "htau_y0_bk_ch", &htau_y0_bk_ch,
 "htau_vh1_bk_ch", &htau_vh1_bk_ch,
 "htau_k1_bk_ch", &htau_k1_bk_ch,
 "htau_vh2_bk_ch", &htau_vh2_bk_ch,
 "htau_k2_bk_ch", &htau_k2_bk_ch,
 "minf_bk_ch", &minf_bk_ch,
 "mtau_bk_ch", &mtau_bk_ch,
 "hinf_bk_ch", &hinf_bk_ch,
 "htau_bk_ch", &htau_bk_ch,
 "zinf_bk_ch", &zinf_bk_ch,
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
 
#define _cvode_ieq _ppvar[4]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"bk_ch",
 "gbar_bk_ch",
 0,
 "g_bk_ch",
 "ik_bk_ch",
 0,
 "m_bk_ch",
 "z_bk_ch",
 "h_bk_ch",
 0,
 0};
 static Symbol* _k_sym;
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 12, _prop);
 	/*initialize range parameters*/
 	gbar = 0.007;
 	_prop->param = _p;
 	_prop->param_size = 12;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[3]._pval = &prop_ion->param[1]; /* cai */
 
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

 void _bk_ch_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", -10000.);
 	ion_reg("ca", -10000.);
 	_k_sym = hoc_lookup("k_ion");
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 12, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 bk_ch C:/Users/bo.bekkouche/PycharmProjects/currentinjection/Snudda/examples/bgd01/parkinson/20211105/PD0/mechanisms/bk_ch.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[3], _dlist1[3];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dz = ( zinf - z ) / ztau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dz = Dz  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ztau )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0 ) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0 ) ) ) / htau ) - h) ;
    z = z + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ztau)))*(- ( ( ( zinf ) ) / ztau ) / ( ( ( ( - 1.0 ) ) ) / ztau ) - z) ;
   }
  return 0;
}
 
static int  rates (  double _lVm ) {
   double _lv ;
 _lv = _lVm + 5.0 ;
   minf = 1.0 / ( 1.0 + exp ( - ( _lv - ( m_vh ) ) / m_k ) ) ;
   mtau = ( 1e3 ) * ( mtau_y0 + 1.0 / ( exp ( ( _lv + mtau_vh1 ) / mtau_k1 ) + exp ( ( _lv + mtau_vh2 ) / mtau_k2 ) ) ) ;
   zinf = 1.0 / ( 1.0 + z_coef / cai ) ;
   hinf = h_y0 + ( 1.0 - h_y0 ) / ( 1.0 + exp ( ( _lv - h_vh ) / h_k ) ) ;
   htau = ( 1e3 ) * ( htau_y0 + 1.0 / ( exp ( ( _lv + htau_vh1 ) / htau_k1 ) + exp ( ( _lv + htau_vh2 ) / htau_k2 ) ) ) ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 3;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
     _ode_spec1 ();
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 3; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ek = _ion_ek;
  cai = _ion_cai;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 1);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  z = z0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   z = zinf ;
   h = hinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  ek = _ion_ek;
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gbar * m * m * m * z * z * h ;
   ik = g * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  cai = _ion_cai;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
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
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
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
  cai = _ion_cai;
 { error =  states();
 if(error){fprintf(stderr,"at line 66 in file bk_ch.mod:\n     SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(z) - _p;  _dlist1[2] = &(Dz) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "bk_ch.mod";
static const char* nmodl_file_text = 
  ": BK-type Purkinje calcium-activated potassium current\n"
  ": Created 8/19/02 - nwg\n"
  "\n"
  "NEURON {\n"
  "       SUFFIX bk_ch\n"
  "       USEION k READ ek WRITE ik\n"
  "       USEION ca READ cai\n"
  "       RANGE gbar, ik, g\n"
  "       GLOBAL minf, mtau, hinf, htau, zinf, ztau\n"
  "       GLOBAL m_vh, m_k, mtau_y0, mtau_vh1, mtau_vh2, mtau_k1, mtau_k2\n"
  "       GLOBAL z_coef, ztau\n"
  "       GLOBAL h_y0, h_vh, h_k, htau_y0, htau_vh1, htau_vh2, htau_k1, htau_k2\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "      (mV) = (millivolt)\n"
  "      (mA) = (milliamp)\n"
  "      (mM) = (milli/liter)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "          v            (mV)\n"
  "          gbar = .007 (mho/cm2)\n"
  "		      \n"
  "          m_vh = -28.9           (mV)\n"
  "          m_k = 6.2            (mV)\n"
  "          mtau_y0 = .000505     (s)\n"
  "          mtau_vh1 = -33.3     (mV)\n"
  "          mtau_k1 = -10         (mV)\n"
  "          mtau_vh2 = 86.4       (mV)\n"
  "          mtau_k2 = 10.1        (mV)\n"
  "\n"
  "          z_coef = .001        (mM)\n"
  "          ztau = 1              (ms)\n"
  "\n"
  "          h_y0 = .085\n"
  "          h_vh = -32          (mV)\n"
  "          h_k = 5.8             (mV)\n"
  "          htau_y0 = .0019      (s)\n"
  "          htau_vh1 = -54.2       (mV)\n"
  "          htau_k1 = -12.9       (mV)\n"
  "          htau_vh2 = 48.5      (mV)\n"
  "          htau_k2 = 5.2        (mV)\n"
  "\n"
  "          ek           (mV)\n"
  "          cai          (mM)\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "         minf\n"
  "         mtau          (ms)\n"
  "         hinf\n"
  "         htau          (ms)\n"
  "         zinf\n"
  "	 g	(S/cm2)\n"
  "         ik            (mA/cm2)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "      m   FROM 0 TO 1\n"
  "      z   FROM 0 TO 1\n"
  "      h   FROM 0 TO 1\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "     SOLVE states METHOD cnexp\n"
  "	   g = gbar * m * m * m * z * z * h\n"
  "           ik = g * (v - ek)\n"
  "}\n"
  "\n"
  "DERIVATIVE states {\n"
  "        rates(v)\n"
  "\n"
  "        m' = (minf - m) / mtau\n"
  "        h' = (hinf - h) / htau\n"
  "        z' = (zinf - z) / ztau\n"
  "}\n"
  "\n"
  "PROCEDURE rates(Vm (mV)) {\n"
  "          LOCAL v\n"
  "          v = Vm + 5\n"
  "          minf = 1 / (1 + exp(-(v - (m_vh)) / m_k))\n"
  "          mtau = (1e3) * (mtau_y0 + 1/(exp((v+ mtau_vh1)/mtau_k1) + exp((v+mtau_vh2)/mtau_k2)))\n"
  "          zinf = 1/(1 + z_coef / cai)\n"
  "          hinf = h_y0 + (1-h_y0) / (1+exp((v - h_vh)/h_k))\n"
  "          htau = (1e3) * (htau_y0 + 1/(exp((v + htau_vh1)/htau_k1)+exp((v+htau_vh2)/htau_k2)))\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "        rates(v)\n"
  "        m = minf\n"
  "        z = zinf\n"
  "        h = hinf\n"
  "}\n"
  "\n"
  "COMMENT\n"
  "From model:\n"
  "Maurice N, Mercer J, Chan CS, Hernandez-Lopez S, Held J, Tkatch T, Surmeier DJ. \n"
  "D2 dopamine receptor-mediated modulation of voltage-dependent Na+ channels \n"
  "reduces autonomous activity in striatal cholinergic interneurons. \n"
  "J Neurosci. 2004 Nov 17;24(46):10289-301. doi: 10.1523/JNEUROSCI.2155-04.2004. \n"
  "PMID: 15548642; PMCID: PMC6730305.\n"
  "\n"
  "\n"
  "ENDCOMMENT\n"
  ;
#endif
