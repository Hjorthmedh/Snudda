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
 
#define nrn_init _nrn_init__hcn12_ch
#define _nrn_initial _nrn_initial__hcn12_ch
#define nrn_cur _nrn_cur__hcn12_ch
#define _nrn_current _nrn_current__hcn12_ch
#define nrn_jacob _nrn_jacob__hcn12_ch
#define nrn_state _nrn_state__hcn12_ch
#define _net_receive _net_receive__hcn12_ch 
#define kin kin__hcn12_ch 
#define rates rates__hcn12_ch 
 
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
#define gbar _p[0]
#define ehcn _p[1]
#define damod _p[2]
#define maxMod _p[3]
#define g _p[4]
#define i _p[5]
#define c _p[6]
#define cac _p[7]
#define o _p[8]
#define cao _p[9]
#define alpha _p[10]
#define beta _p[11]
#define alphaa _p[12]
#define betaa _p[13]
#define Dc _p[14]
#define Dcac _p[15]
#define Do _p[16]
#define Dcao _p[17]
#define v _p[18]
#define _g _p[19]
 
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
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_modulation(void);
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
 "setdata_hcn12_ch", _hoc_setdata,
 "modulation_hcn12_ch", _hoc_modulation,
 "rates_hcn12_ch", _hoc_rates,
 0, 0
};
#define modulation modulation_hcn12_ch
 extern double modulation( _threadargsproto_ );
 /* declare global and static user variables */
#define ai ai_hcn12_ch
 double ai = 1e-05;
#define aac aac_hcn12_ch
 double aac = -0.075;
#define aah aah_hcn12_ch
 double aah = -94.2;
#define aa0 aa0_hcn12_ch
 double aa0 = 0.0006;
#define ac ac_hcn12_ch
 double ac = -0.155;
#define ah ah_hcn12_ch
 double ah = -96;
#define a0 a0_hcn12_ch
 double a0 = 0.006;
#define bf bf_hcn12_ch
 double bf = 8.94;
#define b b_hcn12_ch
 double b = 80;
#define bac bac_hcn12_ch
 double bac = 0.144;
#define bah bah_hcn12_ch
 double bah = -35.5;
#define ba0 ba0_hcn12_ch
 double ba0 = 0.004;
#define bc bc_hcn12_ch
 double bc = 0.144;
#define bh bh_hcn12_ch
 double bh = -51.7;
#define b0 b0_hcn12_ch
 double b0 = 0.0008;
#define gca gca_hcn12_ch
 double gca = 1;
#define koff koff_hcn12_ch
 double koff = 4.5e-05;
#define kon kon_hcn12_ch
 double kon = 30;
#define q10a q10a_hcn12_ch
 double q10a = 1.5;
#define q10v q10v_hcn12_ch
 double q10v = 4;
#define shift shift_hcn12_ch
 double shift = -17;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "a0_hcn12_ch", "/ms",
 "b0_hcn12_ch", "/ms",
 "ah_hcn12_ch", "mV",
 "bh_hcn12_ch", "mV",
 "ac_hcn12_ch", "/mV",
 "bc_hcn12_ch", "/mV",
 "aa0_hcn12_ch", "/ms",
 "ba0_hcn12_ch", "/ms",
 "aah_hcn12_ch", "mV",
 "bah_hcn12_ch", "mV",
 "aac_hcn12_ch", "/mV",
 "bac_hcn12_ch", "/mV",
 "kon_hcn12_ch", "/mM-ms",
 "koff_hcn12_ch", "/ms",
 "ai_hcn12_ch", "mM",
 "shift_hcn12_ch", "mV",
 "gbar_hcn12_ch", "S/cm2",
 "ehcn_hcn12_ch", "mV",
 "g_hcn12_ch", "S/cm2",
 "i_hcn12_ch", "mA/cm2",
 0,0
};
 static double cao0 = 0;
 static double cac0 = 0;
 static double c0 = 0;
 static double delta_t = 0.01;
 static double o0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "a0_hcn12_ch", &a0_hcn12_ch,
 "b0_hcn12_ch", &b0_hcn12_ch,
 "ah_hcn12_ch", &ah_hcn12_ch,
 "bh_hcn12_ch", &bh_hcn12_ch,
 "ac_hcn12_ch", &ac_hcn12_ch,
 "bc_hcn12_ch", &bc_hcn12_ch,
 "aa0_hcn12_ch", &aa0_hcn12_ch,
 "ba0_hcn12_ch", &ba0_hcn12_ch,
 "aah_hcn12_ch", &aah_hcn12_ch,
 "bah_hcn12_ch", &bah_hcn12_ch,
 "aac_hcn12_ch", &aac_hcn12_ch,
 "bac_hcn12_ch", &bac_hcn12_ch,
 "kon_hcn12_ch", &kon_hcn12_ch,
 "koff_hcn12_ch", &koff_hcn12_ch,
 "b_hcn12_ch", &b_hcn12_ch,
 "bf_hcn12_ch", &bf_hcn12_ch,
 "ai_hcn12_ch", &ai_hcn12_ch,
 "gca_hcn12_ch", &gca_hcn12_ch,
 "shift_hcn12_ch", &shift_hcn12_ch,
 "q10v_hcn12_ch", &q10v_hcn12_ch,
 "q10a_hcn12_ch", &q10a_hcn12_ch,
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
"hcn12_ch",
 "gbar_hcn12_ch",
 "ehcn_hcn12_ch",
 "damod_hcn12_ch",
 "maxMod_hcn12_ch",
 0,
 "g_hcn12_ch",
 "i_hcn12_ch",
 0,
 "c_hcn12_ch",
 "cac_hcn12_ch",
 "o_hcn12_ch",
 "cao_hcn12_ch",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 20, _prop);
 	/*initialize range parameters*/
 	gbar = 1;
 	ehcn = -20;
 	damod = 0;
 	maxMod = 1;
 	_prop->param = _p;
 	_prop->param_size = 20;
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
 static void _thread_cleanup(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _hcn12_ch_reg() {
	int _vectorized = 1;
  _initlists();
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 3);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 20, 1);
  hoc_register_dparam_semantics(_mechtype, 0, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 hcn12_ch /home/akozlov/doc/hbp-bsp-live-papers-dev-priv/2020/hjorth_et_al_2020/work/Snudda/x86_64/hcn12_ch.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rates(_threadargsprotocomma_ double);
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  static int _cvspth1 = 1;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 0;
 static int _slist1[4], _dlist1[4]; static double *_temp1;
 static int kin();
 
static int kin (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt)
 {int _reset=0;
 {
   double _lqa ;
 double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<4;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 _lqa = pow( q10a , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   rates ( _threadargscomma_ v ) ;
   /* ~ c <-> o ( alpha , beta )*/
 f_flux =  alpha * c ;
 b_flux =  beta * o ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  alpha ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 3 ,2)  -= _term;
 _term =  beta ;
 _MATELM1( 2 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ c <-> cac ( kon * _lqa * ai / bf , koff * _lqa * b / bf )*/
 f_flux =  kon * _lqa * ai / bf * c ;
 b_flux =  koff * _lqa * b / bf * cac ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  kon * _lqa * ai / bf ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  koff * _lqa * b / bf ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ o <-> cao ( kon * _lqa * ai , koff * _lqa )*/
 f_flux =  kon * _lqa * ai * o ;
 b_flux =  koff * _lqa * cao ;
 _RHS1( 3) -= (f_flux - b_flux);
 
 _term =  kon * _lqa * ai ;
 _MATELM1( 3 ,3)  += _term;
 _term =  koff * _lqa ;
 _MATELM1( 3 ,0)  -= _term;
 /*REACTION*/
  /* ~ cac <-> cao ( alphaa , betaa )*/
 f_flux =  alphaa * cac ;
 b_flux =  betaa * cao ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  alphaa ;
 _MATELM1( 1 ,1)  += _term;
 _term =  betaa ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
   /* c + cac + o + cao = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= cao ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= o ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= cac ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= c ;
 /*CONSERVATION*/
   } return _reset;
 }
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   double _lqv ;
 _lqv = pow( q10v , ( ( celsius - 22.0 ) / 10.0 ) ) ;
   if ( _lv > - 200.0 ) {
     alpha = a0 * _lqv / ( 1.0 + exp ( - ( _lv - ah - shift ) * ac ) ) ;
     beta = b0 * _lqv / ( 1.0 + exp ( - ( _lv - bh - shift ) * bc ) ) ;
     alphaa = aa0 * _lqv / ( 1.0 + exp ( - ( _lv - aah - shift ) * aac ) ) ;
     betaa = ba0 * _lqv / ( 1.0 + exp ( - ( _lv - bah - shift ) * bac ) ) ;
     }
   else {
     alpha = a0 * _lqv / ( 1.0 + exp ( - ( ( - 200.0 ) - ah - shift ) * ac ) ) ;
     beta = b0 * _lqv / ( 1.0 + exp ( - ( ( - 200.0 ) - bh - shift ) * bc ) ) ;
     alphaa = aa0 * _lqv / ( 1.0 + exp ( - ( ( - 200.0 ) - aah - shift ) * aac ) ) ;
     betaa = ba0 * _lqv / ( 1.0 + exp ( - ( ( - 200.0 ) - bah - shift ) * bac ) ) ;
     }
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double modulation ( _threadargsproto_ ) {
   double _lmodulation;
 _lmodulation = 1.0 + damod * ( maxMod - 1.0 ) ;
   
return _lmodulation;
 }
 
static void _hoc_modulation(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  modulation ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double _lqa ;
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<4;_i++) _p[_dlist1[_i]] = 0.0;}
 _lqa = pow( q10a , ( ( celsius - 22.0 ) / 10.0 ) ) ;
 rates ( _threadargscomma_ v ) ;
 /* ~ c <-> o ( alpha , beta )*/
 f_flux =  alpha * c ;
 b_flux =  beta * o ;
 Dc -= (f_flux - b_flux);
 Do += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c <-> cac ( kon * _lqa * ai / bf , koff * _lqa * b / bf )*/
 f_flux =  kon * _lqa * ai / bf * c ;
 b_flux =  koff * _lqa * b / bf * cac ;
 Dc -= (f_flux - b_flux);
 Dcac += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ o <-> cao ( kon * _lqa * ai , koff * _lqa )*/
 f_flux =  kon * _lqa * ai * o ;
 b_flux =  koff * _lqa * cao ;
 Do -= (f_flux - b_flux);
 Dcao += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ cac <-> cao ( alphaa , betaa )*/
 f_flux =  alphaa * cac ;
 b_flux =  betaa * cao ;
 Dcac -= (f_flux - b_flux);
 Dcao += (f_flux - b_flux);
 
 /*REACTION*/
   /* c + cac + o + cao = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double _lqa ;
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<4;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 _lqa = pow( q10a , ( ( celsius - 22.0 ) / 10.0 ) ) ;
 rates ( _threadargscomma_ v ) ;
 /* ~ c <-> o ( alpha , beta )*/
 _term =  alpha ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 3 ,2)  -= _term;
 _term =  beta ;
 _MATELM1( 2 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ c <-> cac ( kon * _lqa * ai / bf , koff * _lqa * b / bf )*/
 _term =  kon * _lqa * ai / bf ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  koff * _lqa * b / bf ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ o <-> cao ( kon * _lqa * ai , koff * _lqa )*/
 _term =  kon * _lqa * ai ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 0 ,3)  -= _term;
 _term =  koff * _lqa ;
 _MATELM1( 3 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ cac <-> cao ( alphaa , betaa )*/
 _term =  alphaa ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  betaa ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
   /* c + cac + o + cao = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 4;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 4; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 4, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  cao = cao0;
  cac = cac0;
  c = c0;
  o = o0;
 {
    _ss_sparse_thread(&_thread[_spth1]._pvoid, 4, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gbar * ( o + cao * gca ) * modulation ( _threadargs_ ) ;
   i = g * ( v - ehcn ) ;
   }
 _current += i;

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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
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
double _dtsav = dt;
if (secondorder) { dt *= 0.5; }
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
 {  sparse_thread(&_thread[_spth1]._pvoid, 4, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 4; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 }}}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(cao) - _p;  _dlist1[0] = &(Dcao) - _p;
 _slist1[1] = &(cac) - _p;  _dlist1[1] = &(Dcac) - _p;
 _slist1[2] = &(c) - _p;  _dlist1[2] = &(Dc) - _p;
 _slist1[3] = &(o) - _p;  _dlist1[3] = &(Do) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/home/akozlov/doc/hbp-bsp-live-papers-dev-priv/2020/hjorth_et_al_2020/work/Snudda/snudda/data/cellspecs/mechanisms/hcn12_ch.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "\n"
  "Josh Held's adaptation to suit HCN1+2.  12/22/2003\n"
  "\n"
  "****\n"
  "Kinetic model of HCN2 channel gating from Wang et al 2002.\n"
  "\n"
  "In this model channel opening is coupled to a change in the affinity of the cyclic nucleotide binding domain for cAMP which is manifest as a shift in the activation curve toward more positive potentials.  This model explains the slow activation kinetics of Ih associated with low concentrations of cAMP.\n"
  "\n"
  "For further details email Matt Nolan at mfnolan@fido.cpmc.columbia.edu.\n"
  "\n"
  "Reference\n"
  "\n"
  "Wang J., Chen S., Nolan M.F. and Siegelbaum S.A. (2002). Activity-dependent regulation of HCN pacemaker channels by cyclicAMP: signalling through dynamic allosteric coupling. Neuron 36, 1-20.\n"
  "****\n"
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
  "	SUFFIX hcn12_ch\n"
  "	NONSPECIFIC_CURRENT i\n"
  "	RANGE i, ehcn, g, gbar\n"
  "	GLOBAL a0, b0, ah, bh, ac, bc, aa0, ba0\n"
  "	GLOBAL aa0, ba0, aah, bah, aac, bac\n"
  "	GLOBAL kon, koff, b, bf, ai, gca, shift\n"
  "    RANGE damod, maxMod\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	(mA) = (milliamp)\n"
  "	(S) = (siemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gbar    = 1		(S/cm2)\n"
  "	ehcn    = -20		(mV)\n"
  "	a0      = .006		(/ms)	: parameters for alpha and beta\n"
  "	b0      = .0008		(/ms)\n"
  "	ah      = -96		(mV)\n"
  "	bh      = -51.7		(mV)\n"
  "	ac      = -.155		(/mV)\n"
  "	bc      = .144		(/mV)\n"
  "	aa0     = .0006		(/ms)	: parameters for alphaa and betaa\n"
  "	ba0     = .004		(/ms)\n"
  "	aah     = -94.2		(mV)\n"
  "	bah     = -35.5		(mV)\n"
  "	aac     = -.075		(/mV)\n"
  "	bac     = .144		(/mV)\n"
  "	kon     = 30		(/mM-ms) : cyclic AMP binding parameters\n"
  "	koff    = 4.5e-05	(/ms)\n"
  "	b       = 80\n"
  "	bf      = 8.94\n"
  "	ai	= 1e-05		(mM)	:concentration cyclic AMP\n"
  "	gca     = 1			: relative conductance of the bound state\n"
  "	shift   = -17		(mV)	: shift in voltage dependence\n"
  "	q10v    = 4                     : q10 value from Magee 1998\n"
  "	q10a    = 1.5			: estimated q10 for the cAMP binding reaction\n"
  "	celsius			(degC)\n"
  "    damod = 0\n"
  "    maxMod = 1\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV)\n"
  "	g	(S/cm2)\n"
  "	i	(mA/cm2)\n"
  "	alpha	(/ms)\n"
  "	beta    (/ms)\n"
  "	alphaa	(/ms)\n"
  "	betaa	(/ms)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	c\n"
  "	cac\n"
  "	o\n"
  "	cao\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	SOLVE kin STEADYSTATE sparse\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE kin METHOD sparse\n"
  "	g = gbar*(o + cao*gca)*modulation()\n"
  "	i = g*(v-ehcn)\n"
  "}\n"
  "\n"
  "KINETIC kin {\n"
  "	LOCAL qa\n"
  "	qa = q10a^((celsius-22 (degC))/10 (degC))\n"
  "	rates(v)\n"
  "	~ c <-> o       (alpha, beta)\n"
  "	~ c <-> cac     (kon*qa*ai/bf,koff*qa*b/bf)\n"
  "	~ o <-> cao     (kon*qa*ai, koff*qa)\n"
  "	~ cac <-> cao   (alphaa, betaa)\n"
  "	CONSERVE c + cac + o + cao = 1\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v(mV)) {\n"
  "	LOCAL qv\n"
  "	qv = q10v^((celsius-22 (degC))/10 (degC))\n"
  "	if (v > -200) {\n"
  "		alpha = a0*qv / (1 + exp(-(v-ah-shift)*ac))\n"
  "		beta = b0*qv / (1 + exp(-(v-bh-shift)*bc))\n"
  "		alphaa = aa0*qv / (1 + exp(-(v-aah-shift)*aac))\n"
  "		betaa = ba0*qv / (1 + exp(-(v-bah-shift)*bac))\n"
  "	} else {\n"
  "		alpha = a0*qv / (1 + exp(-((-200)-ah-shift)*ac))\n"
  "		beta = b0*qv / (1 + exp(-((-200)-bh-shift)*bc))\n"
  "		alphaa = aa0*qv / (1 + exp(-((-200)-aah-shift)*aac))\n"
  "		betaa = ba0*qv / (1 + exp(-((-200)-bah-shift)*bac))\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION modulation() {\n"
  "    : returns modulation factor\n"
  "    \n"
  "    modulation = 1 + damod*(maxMod-1)\n"
  "}\n"
  ;
#endif
