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
 
#define nrn_init _nrn_init__cal_ch
#define _nrn_initial _nrn_initial__cal_ch
#define nrn_cur _nrn_cur__cal_ch
#define _nrn_current _nrn_current__cal_ch
#define nrn_jacob _nrn_jacob__cal_ch
#define nrn_state _nrn_state__cal_ch
#define _net_receive _net_receive__cal_ch 
#define rate rate__cal_ch 
#define state state__cal_ch 
 
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
#define modACh _p[1]
#define maxModACh _p[2]
#define levelACh _p[3]
#define ica _p[4]
#define gcal _p[5]
#define m _p[6]
#define cai _p[7]
#define cao _p[8]
#define Dm _p[9]
#define _g _p[10]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_cao	*_ppvar[1]._pval
#define _ion_ica	*_ppvar[2]._pval
#define _ion_dicadv	*_ppvar[3]._pval
 
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
 static void _hoc_KTF(void);
 static void _hoc_alp(void);
 static void _hoc_bet(void);
 static void _hoc_efun(void);
 static void _hoc_ghk(void);
 static void _hoc_h2(void);
 static void _hoc_modulationACh(void);
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
 "setdata_cal_ch", _hoc_setdata,
 "KTF_cal_ch", _hoc_KTF,
 "alp_cal_ch", _hoc_alp,
 "bet_cal_ch", _hoc_bet,
 "efun_cal_ch", _hoc_efun,
 "ghk_cal_ch", _hoc_ghk,
 "h2_cal_ch", _hoc_h2,
 "modulationACh_cal_ch", _hoc_modulationACh,
 "rate_cal_ch", _hoc_rate,
 0, 0
};
#define KTF KTF_cal_ch
#define _f_bet _f_bet_cal_ch
#define _f_alp _f_alp_cal_ch
#define alp alp_cal_ch
#define bet bet_cal_ch
#define efun efun_cal_ch
#define ghk ghk_cal_ch
#define h2 h2_cal_ch
#define modulationACh modulationACh_cal_ch
 extern double KTF( double );
 extern double _f_bet( double );
 extern double _f_alp( double );
 extern double alp( double );
 extern double bet( double );
 extern double efun( double );
 extern double ghk( double , double , double );
 extern double h2( double );
 extern double modulationACh( );
 /* declare global and static user variables */
#define Ctm Ctm_cal_ch
 double Ctm = 3;
#define atm atm_cal_ch
 double atm = 12;
#define btm btm_cal_ch
 double btm = 11;
#define ki ki_cal_ch
 double ki = 0.001;
#define minf minf_cal_ch
 double minf = 0;
#define tfa tfa_cal_ch
 double tfa = 1;
#define tau tau_cal_ch
 double tau = 0;
#define tm0 tm0_cal_ch
 double tm0 = 0;
#define usetable usetable_cal_ch
 double usetable = 1;
#define vhtm vhtm_cal_ch
 double vhtm = -2;
#define vcm vcm_cal_ch
 double vcm = 5.6;
#define vhm vhm_cal_ch
 double vhm = -1.5;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_cal_ch", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ki_cal_ch", "mM",
 "tau_cal_ch", "ms",
 "gbar_cal_ch", "mho/cm2",
 "ica_cal_ch", "mA/cm2",
 "gcal_cal_ch", "mho/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double m0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ki_cal_ch", &ki_cal_ch,
 "tfa_cal_ch", &tfa_cal_ch,
 "vhm_cal_ch", &vhm_cal_ch,
 "vcm_cal_ch", &vcm_cal_ch,
 "Ctm_cal_ch", &Ctm_cal_ch,
 "atm_cal_ch", &atm_cal_ch,
 "btm_cal_ch", &btm_cal_ch,
 "tm0_cal_ch", &tm0_cal_ch,
 "vhtm_cal_ch", &vhtm_cal_ch,
 "minf_cal_ch", &minf_cal_ch,
 "tau_cal_ch", &tau_cal_ch,
 "usetable_cal_ch", &usetable_cal_ch,
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
"cal_ch",
 "gbar_cal_ch",
 "modACh_cal_ch",
 "maxModACh_cal_ch",
 "levelACh_cal_ch",
 0,
 "ica_cal_ch",
 "gcal_cal_ch",
 0,
 "m_cal_ch",
 0,
 0};
 static Symbol* _ca_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	gbar = 0.003;
 	modACh = 0;
 	maxModACh = 1;
 	levelACh = 0;
 	_prop->param = _p;
 	_prop->param_size = 11;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 5, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 	_ppvar[1]._pval = &prop_ion->param[2]; /* cao */
 	_ppvar[2]._pval = &prop_ion->param[3]; /* ica */
 	_ppvar[3]._pval = &prop_ion->param[4]; /* _ion_dicadv */
 
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

 void _cal_ch_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 11, 5);
  hoc_register_dparam_semantics(_mechtype, 0, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "ca_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 cal_ch C:/Users/bo.bekkouche/PycharmProjects/currentinjection/Snudda/examples/bgd01/parkinson/20211105/PD0/mechanisms/cal_ch.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 
#define FARADAY _nrnunit_FARADAY[_nrnunit_use_legacy_]
static double _nrnunit_FARADAY[2] = {0xc.0f87d73baa128p+3, 96.4853}; /* 96.4853321233100161 */
 
#define R _nrnunit_R[_nrnunit_use_legacy_]
static double _nrnunit_R[2] = {0x8.50809f44c85fp+0, 8.3145}; /* 8.3144626181532395 */
 static double KTOMV = .0853;
 static double *_t_alp;
 static double *_t_bet;
static int _reset;
static char *modelname = "l-calcium channel";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int rate(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[1], _dlist1[1];
 static int state(_threadargsproto_);
 static double _n_bet(double);
 static double _n_alp(double);
 
double h2 (  double _lcai ) {
   double _lh2;
 _lh2 = ki / ( ki + _lcai ) ;
   
return _lh2;
 }
 
static void _hoc_h2(void) {
  double _r;
   _r =  h2 (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double ghk (  double _lv , double _lci , double _lco ) {
   double _lghk;
 double _lnu , _lf ;
 _lf = KTF ( _threadargscomma_ celsius ) / 2.0 ;
   _lnu = _lv / _lf ;
   _lghk = - _lf * ( 1. - ( _lci / _lco ) * exp ( _lnu ) ) * efun ( _threadargscomma_ _lnu ) ;
   
return _lghk;
 }
 
static void _hoc_ghk(void) {
  double _r;
   _r =  ghk (  *getarg(1) , *getarg(2) , *getarg(3) );
 hoc_retpushx(_r);
}
 
double KTF (  double _lcelsius ) {
   double _lKTF;
 _lKTF = ( ( 25. / 293.15 ) * ( _lcelsius + 273.15 ) ) ;
   
return _lKTF;
 }
 
static void _hoc_KTF(void) {
  double _r;
   _r =  KTF (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double efun (  double _lz ) {
   double _lefun;
 if ( fabs ( _lz ) < 1e-4 ) {
     _lefun = 1.0 - _lz / 2.0 ;
     }
   else {
     _lefun = _lz / ( exp ( _lz ) - 1.0 ) ;
     }
   
return _lefun;
 }
 
static void _hoc_efun(void) {
  double _r;
   _r =  efun (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_alp, _tmin_alp;
 static void _check_alp();
 static void _check_alp() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_alp =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_alp)/200.; _mfac_alp = 1./_dx;
   for (_i=0, _x=_tmin_alp; _i < 201; _x += _dx, _i++) {
    _t_alp[_i] = _f_alp(_x);
   }
  }
 }

 double alp(double _lv){ _check_alp();
 return _n_alp(_lv);
 }

 static double _n_alp(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_alp(_lv); 
}
 _xi = _mfac_alp * (_lv - _tmin_alp);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_alp[0];
 }
 if (_xi >= 200.) {
 return _t_alp[200];
 }
 _i = (int) _xi;
 return _t_alp[_i] + (_xi - (double)_i)*(_t_alp[_i+1] - _t_alp[_i]);
 }

 
double _f_alp (  double _lv ) {
   double _lalp;
 _lalp = 15.69 * ( - 1.0 * _lv + 81.5 ) / ( exp ( ( - 1.0 * _lv + 81.5 ) / 10.0 ) - 1.0 ) ;
   
return _lalp;
 }
 
static void _hoc_alp(void) {
  double _r;
    _r =  alp (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_bet, _tmin_bet;
 static void _check_bet();
 static void _check_bet() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  if (!usetable) {return;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_bet =  - 150.0 ;
   _tmax =  150.0 ;
   _dx = (_tmax - _tmin_bet)/200.; _mfac_bet = 1./_dx;
   for (_i=0, _x=_tmin_bet; _i < 201; _x += _dx, _i++) {
    _t_bet[_i] = _f_bet(_x);
   }
  }
 }

 double bet(double _lv){ _check_bet();
 return _n_bet(_lv);
 }

 static double _n_bet(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 return _f_bet(_lv); 
}
 _xi = _mfac_bet * (_lv - _tmin_bet);
 if (isnan(_xi)) {
  return _xi; }
 if (_xi <= 0.) {
 return _t_bet[0];
 }
 if (_xi >= 200.) {
 return _t_bet[200];
 }
 _i = (int) _xi;
 return _t_bet[_i] + (_xi - (double)_i)*(_t_bet[_i+1] - _t_bet[_i]);
 }

 
double _f_bet (  double _lv ) {
   double _lbet;
 _lbet = 0.29 * exp ( - _lv / 10.86 ) ;
   
return _lbet;
 }
 
static void _hoc_bet(void) {
  double _r;
    _r =  bet (  *getarg(1) );
 hoc_retpushx(_r);
}
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   rate ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / tau ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 rate ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau )) ;
  return 0;
}
 /*END CVODE*/
 static int state () {_reset=0;
 {
   rate ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau)))*(- ( ( ( minf ) ) / tau ) / ( ( ( ( - 1.0 ) ) ) / tau ) - m) ;
   }
  return 0;
}
 
static int  rate (  double _lv ) {
   double _la ;
 _la = alp ( _threadargscomma_ _lv ) ;
   tau = Ctm / ( exp ( ( _lv - vhtm ) / atm ) + exp ( ( vhtm - _lv ) / btm ) ) + tm0 ;
   minf = 1.0 / ( 1.0 + exp ( - ( _lv - vhm ) / vcm ) ) ;
    return 0; }
 
static void _hoc_rate(void) {
  double _r;
   _r = 1.;
 rate (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double modulationACh (  ) {
   double _lmodulationACh;
 _lmodulationACh = 1.0 + modACh * ( maxModACh - 1.0 ) * levelACh ;
   
return _lmodulationACh;
 }
 
static void _hoc_modulationACh(void) {
  double _r;
   _r =  modulationACh (  );
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
  cai = _ion_cai;
  cao = _ion_cao;
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
  cai = _ion_cai;
  cao = _ion_cao;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 1, 2);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 2, 3);
   nrn_update_ion_pointer(_ca_sym, _ppvar, 3, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  m = m0;
 {
   rate ( _threadargscomma_ v ) ;
   m = minf ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   gcal = gbar * m * m * h2 ( _threadargscomma_ cai ) * modulationACh ( _threadargs_ ) ;
   ica = gcal * ghk ( _threadargscomma_ v , cai , cao ) ;
   }
 _current += ica;

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
  cai = _ion_cai;
  cao = _ion_cao;
 _g = _nrn_current(_v + .001);
 	{ double _dica;
  _dica = ica;
 _rhs = _nrn_current(_v);
  _ion_dicadv += (_dica - ica)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ica += ica ;
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
  cai = _ion_cai;
  cao = _ion_cao;
 { error =  state();
 if(error){fprintf(stderr,"at line 86 in file cal_ch.mod:\n	SOLVE state METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_alp = makevector(201*sizeof(double));
   _t_bet = makevector(201*sizeof(double));
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "cal_ch.mod";
static const char* nmodl_file_text = 
  ":Migliore file Modify by Maciej Lazarewicz (mailto:mlazarew@gmu.edu) May/16/2001\n"
  "\n"
  "TITLE l-calcium channel\n"
  ": l-type calcium channel\n"
  "\n"
  ": copy by josh\n"
  "\n"
  "COMMENT\n"
  "\n"
  "Neuromodulation is added as functions:\n"
  "    \n"
  "    modulationDA = 1 + modDA*(maxModDA-1)*levelDA\n"
  "\n"
  "where:\n"
  "    \n"
  "    modDA  [0]: is a switch for turning modulation on or off {1/0}\n"
  "    maxModDA [1]: is the maximum modulation for this specific channel (read from the param file)\n"
  "                    e.g. 10% increase would correspond to a factor of 1.1 (100% +10%) {0-inf}\n"
  "    levelDA  [0]: is an additional parameter for scaling modulation. \n"
  "                Can be used simulate non static modulation by gradually changing the value from 0 to 1 {0-1}\n"
  "									\n"
  "	  Further neuromodulators can be added by for example:\n"
  "          modulationDA = 1 + modDA*(maxModDA-1)\n"
  "	  modulationACh = 1 + modACh*(maxModACh-1)\n"
  "	  ....\n"
  "\n"
  "	  etc. for other neuromodulators\n"
  "	  \n"
  "	   \n"
  "								     \n"
  "[] == default values\n"
  "{} == ranges\n"
  "    \n"
  "ENDCOMMENT\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX cal_ch\n"
  "	USEION ca READ cai,cao WRITE ica\n"
  "        RANGE  gbar,ica , gcal\n"
  "	GLOBAL vhm, vcm\n"
  "	GLOBAL Ctm, atm, btm, tm0, vhtm\n"
  "        GLOBAL minf,tau\n"
  "        RANGE modACh, maxModACh, levelACh\n"
  "\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mA) 	= 	(milliamp)\n"
  "	(mV) 	= 	(millivolt)\n"
  "	FARADAY =  	(faraday)  (kilocoulombs)\n"
  "	R 	= 	(k-mole) (joule/degC)\n"
  "	KTOMV 	= .0853 (mV/degC)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	v (mV)\n"
  "	celsius = 6.3	(degC)\n"
  "	gbar	= .003 	(mho/cm2)\n"
  "	ki	= .001 	(mM)\n"
  "	cai 		(mM)\n"
  "	cao 		(mM)\n"
  "        tfa	= 1\n"
  "	vhm = -1.5\n"
  "	vcm = 5.6\n"
  "	Ctm = 3\n"
  "	atm = 12\n"
  "	btm = 11\n"
  "	tm0 = 0\n"
  "	vhtm = -2\n"
  "        modACh = 0\n"
  "        maxModACh = 1\n"
  "        levelACh = 0\n"
  "\n"
  "}\n"
  "\n"
  "STATE { m }\n"
  "\n"
  "ASSIGNED {\n"
  "	ica (mA/cm2)\n"
  "        gcal (mho/cm2)\n"
  "        minf\n"
  "        tau   (ms)\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	gcal = gbar*m*m*h2(cai)*modulationACh()\n"
  "	ica  = gcal*ghk(v,cai,cao)\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	rate(v)\n"
  "	m = minf\n"
  "}\n"
  "\n"
  "FUNCTION h2(cai(mM)) {\n"
  "	h2 = ki/(ki+cai)\n"
  "}\n"
  "\n"
  "FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {\n"
  "        LOCAL nu,f\n"
  "\n"
  "        f = KTF(celsius)/2\n"
  "        nu = v/f\n"
  "        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)\n"
  "}\n"
  "\n"
  "FUNCTION KTF(celsius (DegC)) (mV) {\n"
  "        KTF = ((25./293.15)*(celsius + 273.15))\n"
  "}\n"
  "\n"
  "FUNCTION efun(z) {\n"
  "	if (fabs(z) < 1e-4) {\n"
  "		efun = 1 - z/2\n"
  "	}else{\n"
  "		efun = z/(exp(z) - 1)\n"
  "	}\n"
  "}\n"
  "\n"
  "FUNCTION alp(v(mV)) (1/ms) {\n"
  "	TABLE FROM -150 TO 150 WITH 200\n"
  "	alp = 15.69*(-1.0*v+81.5)/(exp((-1.0*v+81.5)/10.0)-1.0)\n"
  "}\n"
  "\n"
  "FUNCTION bet(v(mV)) (1/ms) {\n"
  "	TABLE FROM -150 TO 150 WITH 200\n"
  "	bet = 0.29*exp(-v/10.86)\n"
  "}\n"
  "\n"
  "DERIVATIVE state {  \n"
  "        rate(v)\n"
  "        m' = (minf - m)/tau\n"
  "}\n"
  "\n"
  "PROCEDURE rate(v (mV)) { :callable from hoc\n"
  "        LOCAL a\n"
  "\n"
  "        a    = alp(v)\n"
  "        :tau  = 1/(tfa*(a + bet(v)))\n"
  "        :minf = tfa*a*tau\n"
  "	tau = Ctm/(exp((v-vhtm)/atm) + exp((vhtm-v)/btm)) + tm0\n"
  "	minf = 1/(1+exp(-(v-vhm)/vcm))\n"
  "}\n"
  " \n"
  "FUNCTION modulationACh() {\n"
  "    : returns modulation factor\n"
  "    \n"
  "    modulationACh = 1 + modACh*(maxModACh-1)*levelACh \n"
  "}\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  "\n"
  ;
#endif
