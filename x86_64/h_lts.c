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
 
#define nrn_init _nrn_init__hd_lts
#define _nrn_initial _nrn_initial__hd_lts
#define nrn_cur _nrn_cur__hd_lts
#define _nrn_current _nrn_current__hd_lts
#define nrn_jacob _nrn_jacob__hd_lts
#define nrn_state _nrn_state__hd_lts
#define _net_receive _net_receive__hd_lts 
#define rate rate__hd_lts 
#define states states__hd_lts 
 
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
#define ghdbar _p[0]
#define vhalfl _p[1]
#define damod _p[2]
#define maxMod _p[3]
#define i _p[4]
#define l _p[5]
#define Dl _p[6]
#define ghd _p[7]
#define _g _p[8]
 
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
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_alpt(void);
 static void _hoc_bett(void);
 static void _hoc_modulation(void);
 static void _hoc_rate(void);
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
 "setdata_hd_lts", _hoc_setdata,
 "alpt_hd_lts", _hoc_alpt,
 "bett_hd_lts", _hoc_bett,
 "modulation_hd_lts", _hoc_modulation,
 "rate_hd_lts", _hoc_rate,
 0, 0
};
#define alpt alpt_hd_lts
#define bett bett_hd_lts
#define modulation modulation_hd_lts
 extern double alpt( double );
 extern double bett( double );
 extern double modulation( );
 /* declare global and static user variables */
#define a0t a0t_hd_lts
 double a0t = 0.011;
#define ehd ehd_hd_lts
 double ehd = -30;
#define gmt gmt_hd_lts
 double gmt = 0.4;
#define kl kl_hd_lts
 double kl = -8;
#define linf linf_hd_lts
 double linf = 0;
#define qtl qtl_hd_lts
 double qtl = 1;
#define q10 q10_hd_lts
 double q10 = 4.5;
#define taul taul_hd_lts
 double taul = 0;
#define vhalft vhalft_hd_lts
 double vhalft = -75;
#define zetat zetat_hd_lts
 double zetat = 2.2;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ehd_hd_lts", "mV",
 "vhalft_hd_lts", "mV",
 "a0t_hd_lts", "/ms",
 "zetat_hd_lts", "1",
 "gmt_hd_lts", "1",
 "ghdbar_hd_lts", "mho/cm2",
 "vhalfl_hd_lts", "mV",
 "i_hd_lts", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double l0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ehd_hd_lts", &ehd_hd_lts,
 "kl_hd_lts", &kl_hd_lts,
 "vhalft_hd_lts", &vhalft_hd_lts,
 "a0t_hd_lts", &a0t_hd_lts,
 "zetat_hd_lts", &zetat_hd_lts,
 "gmt_hd_lts", &gmt_hd_lts,
 "q10_hd_lts", &q10_hd_lts,
 "qtl_hd_lts", &qtl_hd_lts,
 "linf_hd_lts", &linf_hd_lts,
 "taul_hd_lts", &taul_hd_lts,
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
 
#define _cvode_ieq _ppvar[0]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"hd_lts",
 "ghdbar_hd_lts",
 "vhalfl_hd_lts",
 "damod_hd_lts",
 "maxMod_hd_lts",
 0,
 "i_hd_lts",
 0,
 "l_hd_lts",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 9, _prop);
 	/*initialize range parameters*/
 	ghdbar = 0.0001;
 	vhalfl = -81;
 	damod = 0;
 	maxMod = 1;
 	_prop->param = _p;
 	_prop->param_size = 9;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 1, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _h_lts_reg() {
	int _vectorized = 0;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 9, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 hd_lts /home/akozlov/doc/hbp-bsp-live-papers-dev-priv/2020/hjorth_et_al_2020/work/Snudda/x86_64/h_lts.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "I-h channel from Magee 1998 for distal dendrites";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int states(_threadargsproto_);
 
double alpt (  double _lv ) {
   double _lalpt;
 _lalpt = exp ( 0.0378 * zetat * ( _lv - vhalft ) ) ;
   
return _lalpt;
 }
 
static void _hoc_alpt(void) {
  double _r;
   _r =  alpt (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double bett (  double _lv ) {
   double _lbett;
 _lbett = exp ( 0.0378 * zetat * gmt * ( _lv - vhalft ) ) ;
   
return _lbett;
 }
 
static void _hoc_bett(void) {
  double _r;
   _r =  bett (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rate ( _threadargscomma_ v ) ;
   Dl = ( linf - l ) / taul ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rate ( _threadargscomma_ v ) ;
 Dl = Dl  / (1. - dt*( ( ( ( - 1.0 ) ) ) / taul )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   rate ( _threadargscomma_ v ) ;
    l = l + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / taul)))*(- ( ( ( linf ) ) / taul ) / ( ( ( ( - 1.0 ) ) ) / taul ) - l) ;
   }
  return 0;
}
 
static int  rate (  double _lv ) {
   double _la , _lqt ;
 _lqt = pow( q10 , ( ( celsius - 33.0 ) / 10.0 ) ) ;
   _la = alpt ( _threadargscomma_ _lv ) ;
   linf = 1.0 / ( 1.0 + exp ( - ( _lv - vhalfl ) / kl ) ) ;
   taul = bett ( _threadargscomma_ _lv ) / ( qtl * _lqt * a0t * ( 1.0 + _la ) ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   _r = 1.;
 rate (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double modulation (  ) {
   double _lmodulation;
 _lmodulation = 1.0 + damod * ( maxMod - 1.0 ) ;
   
return _lmodulation;
 }
 
static void _hoc_modulation(void) {
  double _r;
   _r =  modulation (  );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 1;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 ();
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 1; ++_i) {
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
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  l = l0;
 {
   rate ( _threadargscomma_ v ) ;
   l = linf ;
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
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ghd = ghdbar * l * modulation ( _threadargs_ ) ;
   i = ghd * ( v - ehd ) ;
   }
 _current += i;

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
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
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
 { error =  states();
 if(error){fprintf(stderr,"at line 71 in file h_lts.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }}}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(l) - _p;  _dlist1[0] = &(Dl) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/akozlov/doc/hbp-bsp-live-papers-dev-priv/2020/hjorth_et_al_2020/work/Snudda/snudda/data/cellspecs/mechanisms/h_lts.mod";
static const char* nmodl_file_text = 
  "TITLE I-h channel from Magee 1998 for distal dendrites\n"
  "\n"
  "COMMENT\n"
  "\n"
  "neuromodulation is added as functions:\n"
  "    \n"
  "    modulation = 1 + damod*(maxMod-1)\n"
  "\n"
  "where:\n"
  "    \n"
  "    damod  [0]: is a switch for turning modulation on or off {1/0}\n"
  "    maxMod [1]: is the maximum modulation for this specific channel (read from the param file)\n"
  "                    e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}\n"
  "\n"
  "[] == default values\n"
  "{} == ranges\n"
  "    \n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX hd_lts\n"
  "	NONSPECIFIC_CURRENT i\n"
  "    RANGE ghdbar, vhalfl\n"
  "    GLOBAL linf,taul\n"
  "    RANGE damod, maxMod\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	v 		        (mV)\n"
  "    ehd  = -30		(mV)        \n"
  "	celsius 	    (degC)\n"
  "	ghdbar=.0001 	(mho/cm2)\n"
  "    vhalfl=-81   	(mV)\n"
  "	kl=-8\n"
  "    vhalft=-75   	(mV)\n"
  "    a0t=0.011      	(/ms)\n"
  "    zetat=2.2    	(1)\n"
  "    gmt=.4   	    (1)\n"
  "	q10=4.5\n"
  "	qtl=1\n"
  "    damod = 0\n"
  "    maxMod = 1\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "STATE {\n"
  "    l\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	i (mA/cm2)\n"
  "    linf      \n"
  "    taul\n"
  "    ghd\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rate(v)\n"
  "	l=linf\n"
  "}\n"
  "\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "	ghd = ghdbar*l*modulation()\n"
  "	i = ghd*(v-ehd)\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  "FUNCTION alpt(v(mV)) {\n"
  "    alpt = exp(0.0378*zetat*(v-vhalft)) \n"
  "}\n"
  "\n"
  "FUNCTION bett(v(mV)) {\n"
  "    bett = exp(0.0378*zetat*gmt*(v-vhalft)) \n"
  "}\n"
  "\n"
  "DERIVATIVE states {     : exact when v held constant; integrates over dt step\n"
  "    rate(v)\n"
  "    l' =  (linf - l)/taul\n"
  "}\n"
  "\n"
  "PROCEDURE rate(v (mV)) { :callable from hoc\n"
  "    LOCAL a,qt\n"
  "    qt=q10^((celsius-33)/10)\n"
  "    a = alpt(v)\n"
  "    linf = 1/(1 + exp(-(v-vhalfl)/kl))\n"
  "    :linf = 1/(1+ alpl(v))\n"
  "    taul = bett(v)/(qtl*qt*a0t*(1+a))\n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "FUNCTION modulation() {\n"
  "    : returns modulation factor\n"
  "    \n"
  "    modulation = 1 + damod*(maxMod-1)\n"
  "}\n"
  "\n"
  "\n"
  ;
#endif
