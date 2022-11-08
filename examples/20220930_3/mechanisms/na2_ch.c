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
 
#define nrn_init _nrn_init__na2_ch
#define _nrn_initial _nrn_initial__na2_ch
#define nrn_cur _nrn_cur__na2_ch
#define _nrn_current _nrn_current__na2_ch
#define nrn_jacob _nrn_jacob__na2_ch
#define nrn_state _nrn_state__na2_ch
#define _net_receive _net_receive__na2_ch 
#define kin kin__na2_ch 
#define rates rates__na2_ch 
 
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
#define modDA _p[1]
#define maxModDA _p[2]
#define levelDA _p[3]
#define g _p[4]
#define ina _p[5]
#define a _p[6]
#define c1 _p[7]
#define c2 _p[8]
#define c3 _p[9]
#define c4 _p[10]
#define c5 _p[11]
#define ct _p[12]
#define o _p[13]
#define i1 _p[14]
#define i2 _p[15]
#define i3 _p[16]
#define i4 _p[17]
#define i5 _p[18]
#define i6 _p[19]
#define ift _p[20]
#define is1 _p[21]
#define is2 _p[22]
#define ist _p[23]
#define it _p[24]
#define ena _p[25]
#define alpha _p[26]
#define beta _p[27]
#define gamma _p[28]
#define delta _p[29]
#define Dc1 _p[30]
#define Dc2 _p[31]
#define Dc3 _p[32]
#define Dc4 _p[33]
#define Dc5 _p[34]
#define Dct _p[35]
#define Do _p[36]
#define Di1 _p[37]
#define Di2 _p[38]
#define Di3 _p[39]
#define Di4 _p[40]
#define Di5 _p[41]
#define Di6 _p[42]
#define Dift _p[43]
#define Dis1 _p[44]
#define Dis2 _p[45]
#define Dist _p[46]
#define Dit _p[47]
#define v _p[48]
#define _g _p[49]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
 
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
 static void _hoc_modulationDA(void);
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
 "setdata_na2_ch", _hoc_setdata,
 "modulationDA_na2_ch", _hoc_modulationDA,
 "rates_na2_ch", _hoc_rates,
 0, 0
};
#define modulationDA modulationDA_na2_ch
 extern double modulationDA( _threadargsproto_ );
 /* declare global and static user variables */
#define Coff Coff_na2_ch
 double Coff = 0.1;
#define Con Con_na2_ch
 double Con = 0.001;
#define Ooff Ooff_na2_ch
 double Ooff = 0.01;
#define Oon Oon_na2_ch
 double Oon = 1.6;
#define aS2 aS2_na2_ch
 double aS2 = 0.0002;
#define aS1 aS1_na2_ch
 double aS1 = 0.0025;
#define a0 a0_na2_ch
 double a0 = 37;
#define bS bS_na2_ch
 double bS = 0.00017;
#define b0 b0_na2_ch
 double b0 = 10;
#define d0 d0_na2_ch
 double d0 = 30;
#define g0 g0_na2_ch
 double g0 = 40;
#define vcb vcb_na2_ch
 double vcb = -10;
#define vhb vhb_na2_ch
 double vhb = -50;
#define vca vca_na2_ch
 double vca = 40;
#define vha vha_na2_ch
 double vha = 45;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "a0_na2_ch", "1/ms",
 "vha_na2_ch", "mV",
 "vca_na2_ch", "mV",
 "b0_na2_ch", "1/ms",
 "vhb_na2_ch", "mV",
 "vcb_na2_ch", "mV",
 "g0_na2_ch", "1/ms",
 "d0_na2_ch", "1/ms",
 "aS1_na2_ch", "1/ms",
 "aS2_na2_ch", "1/ms",
 "bS_na2_ch", "1/ms",
 "Con_na2_ch", "1/ms",
 "Coff_na2_ch", "1/ms",
 "Oon_na2_ch", "1/ms",
 "Ooff_na2_ch", "1/ms",
 "gbar_na2_ch", "S/cm2",
 "g_na2_ch", "S/cm2",
 "ina_na2_ch", "mA/cm2",
 0,0
};
 static double ct0 = 0;
 static double c50 = 0;
 static double c40 = 0;
 static double c30 = 0;
 static double c20 = 0;
 static double c10 = 0;
 static double delta_t = 0.01;
 static double it0 = 0;
 static double ist0 = 0;
 static double is20 = 0;
 static double is10 = 0;
 static double ift0 = 0;
 static double i60 = 0;
 static double i50 = 0;
 static double i40 = 0;
 static double i30 = 0;
 static double i20 = 0;
 static double i10 = 0;
 static double o0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "a0_na2_ch", &a0_na2_ch,
 "vha_na2_ch", &vha_na2_ch,
 "vca_na2_ch", &vca_na2_ch,
 "b0_na2_ch", &b0_na2_ch,
 "vhb_na2_ch", &vhb_na2_ch,
 "vcb_na2_ch", &vcb_na2_ch,
 "g0_na2_ch", &g0_na2_ch,
 "d0_na2_ch", &d0_na2_ch,
 "aS1_na2_ch", &aS1_na2_ch,
 "aS2_na2_ch", &aS2_na2_ch,
 "bS_na2_ch", &bS_na2_ch,
 "Con_na2_ch", &Con_na2_ch,
 "Coff_na2_ch", &Coff_na2_ch,
 "Oon_na2_ch", &Oon_na2_ch,
 "Ooff_na2_ch", &Ooff_na2_ch,
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
"na2_ch",
 "gbar_na2_ch",
 "modDA_na2_ch",
 "maxModDA_na2_ch",
 "levelDA_na2_ch",
 0,
 "g_na2_ch",
 "ina_na2_ch",
 "a_na2_ch",
 0,
 "c1_na2_ch",
 "c2_na2_ch",
 "c3_na2_ch",
 "c4_na2_ch",
 "c5_na2_ch",
 "ct_na2_ch",
 "o_na2_ch",
 "i1_na2_ch",
 "i2_na2_ch",
 "i3_na2_ch",
 "i4_na2_ch",
 "i5_na2_ch",
 "i6_na2_ch",
 "ift_na2_ch",
 "is1_na2_ch",
 "is2_na2_ch",
 "ist_na2_ch",
 "it_na2_ch",
 0,
 0};
 static Symbol* _na_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 50, _prop);
 	/*initialize range parameters*/
 	gbar = 1;
 	modDA = 0;
 	maxModDA = 1;
 	levelDA = 0;
 	_prop->param = _p;
 	_prop->param_size = 50;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 4, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _na2_ch_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 3);
  _extcall_thread = (Datum*)ecalloc(2, sizeof(Datum));
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 50, 4);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 na2_ch C:/Users/bo.bekkouche/PycharmProjects/currentinjection/Snudda/examples/bgd01/parkinson/20211105/PD0/mechanisms/na2_ch.mod\n");
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
 extern double *_nrn_thread_getelm();
 
#define _MATELM1(_row,_col) *(_nrn_thread_getelm(_so, _row + 1, _col + 1))
 
#define _RHS1(_arg) _rhs[_arg+1]
  
#define _linmat1  1
 static int _spth1 = 1;
 static int _cvspth1 = 0;
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[14], _dlist1[14]; static double *_temp1;
 static int kin();
 
static int kin (void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt)
 {int _reset=0;
 {
   double b_flux, f_flux, _term; int _i;
 {int _i; double _dt1 = 1.0/dt;
for(_i=1;_i<14;_i++){
  	_RHS1(_i) = -_dt1*(_p[_slist1[_i]] - _p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ v ) ;
   /* ~ c1 <-> c2 ( 4.0 * alpha , beta )*/
 f_flux =  4.0 * alpha * c1 ;
 b_flux =  beta * c2 ;
 _RHS1( 5) -= (f_flux - b_flux);
 _RHS1( 4) += (f_flux - b_flux);
 
 _term =  4.0 * alpha ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 4 ,5)  -= _term;
 _term =  beta ;
 _MATELM1( 5 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ c2 <-> c3 ( 3.0 * alpha , 2.0 * beta )*/
 f_flux =  3.0 * alpha * c2 ;
 b_flux =  2.0 * beta * c3 ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 3) += (f_flux - b_flux);
 
 _term =  3.0 * alpha ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  2.0 * beta ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ c3 <-> c4 ( 2.0 * alpha , 3.0 * beta )*/
 f_flux =  2.0 * alpha * c3 ;
 b_flux =  3.0 * beta * c4 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 2) += (f_flux - b_flux);
 
 _term =  2.0 * alpha ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  3.0 * beta ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ c4 <-> c5 ( alpha , 4.0 * beta )*/
 f_flux =  alpha * c4 ;
 b_flux =  4.0 * beta * c5 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 1) += (f_flux - b_flux);
 
 _term =  alpha ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  4.0 * beta ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ c5 <-> o ( gamma , delta )*/
 f_flux =  gamma * c5 ;
 b_flux =  delta * o ;
 _RHS1( 1) -= (f_flux - b_flux);
 
 _term =  gamma ;
 _MATELM1( 1 ,1)  += _term;
 _term =  delta ;
 _MATELM1( 1 ,0)  -= _term;
 /*REACTION*/
  /* ~ o <-> is1 ( aS1 , bS )*/
 f_flux =  aS1 * o ;
 b_flux =  bS * is1 ;
 _RHS1( 7) += (f_flux - b_flux);
 
 _term =  aS1 ;
 _MATELM1( 7 ,0)  -= _term;
 _term =  bS ;
 _MATELM1( 7 ,7)  += _term;
 /*REACTION*/
  /* ~ i1 <-> i2 ( 4.0 * alpha * a , beta / a )*/
 f_flux =  4.0 * alpha * a * i1 ;
 b_flux =  beta / a * i2 ;
 _RHS1( 13) -= (f_flux - b_flux);
 _RHS1( 12) += (f_flux - b_flux);
 
 _term =  4.0 * alpha * a ;
 _MATELM1( 13 ,13)  += _term;
 _MATELM1( 12 ,13)  -= _term;
 _term =  beta / a ;
 _MATELM1( 13 ,12)  -= _term;
 _MATELM1( 12 ,12)  += _term;
 /*REACTION*/
  /* ~ i2 <-> i3 ( 3.0 * alpha * a , 2.0 * beta / a )*/
 f_flux =  3.0 * alpha * a * i2 ;
 b_flux =  2.0 * beta / a * i3 ;
 _RHS1( 12) -= (f_flux - b_flux);
 _RHS1( 11) += (f_flux - b_flux);
 
 _term =  3.0 * alpha * a ;
 _MATELM1( 12 ,12)  += _term;
 _MATELM1( 11 ,12)  -= _term;
 _term =  2.0 * beta / a ;
 _MATELM1( 12 ,11)  -= _term;
 _MATELM1( 11 ,11)  += _term;
 /*REACTION*/
  /* ~ i3 <-> i4 ( 2.0 * alpha * a , 3.0 * beta / a )*/
 f_flux =  2.0 * alpha * a * i3 ;
 b_flux =  3.0 * beta / a * i4 ;
 _RHS1( 11) -= (f_flux - b_flux);
 _RHS1( 10) += (f_flux - b_flux);
 
 _term =  2.0 * alpha * a ;
 _MATELM1( 11 ,11)  += _term;
 _MATELM1( 10 ,11)  -= _term;
 _term =  3.0 * beta / a ;
 _MATELM1( 11 ,10)  -= _term;
 _MATELM1( 10 ,10)  += _term;
 /*REACTION*/
  /* ~ i4 <-> i5 ( alpha * a , 4.0 * beta / a )*/
 f_flux =  alpha * a * i4 ;
 b_flux =  4.0 * beta / a * i5 ;
 _RHS1( 10) -= (f_flux - b_flux);
 _RHS1( 9) += (f_flux - b_flux);
 
 _term =  alpha * a ;
 _MATELM1( 10 ,10)  += _term;
 _MATELM1( 9 ,10)  -= _term;
 _term =  4.0 * beta / a ;
 _MATELM1( 10 ,9)  -= _term;
 _MATELM1( 9 ,9)  += _term;
 /*REACTION*/
  /* ~ i5 <-> i6 ( gamma , delta )*/
 f_flux =  gamma * i5 ;
 b_flux =  delta * i6 ;
 _RHS1( 9) -= (f_flux - b_flux);
 _RHS1( 8) += (f_flux - b_flux);
 
 _term =  gamma ;
 _MATELM1( 9 ,9)  += _term;
 _MATELM1( 8 ,9)  -= _term;
 _term =  delta ;
 _MATELM1( 9 ,8)  -= _term;
 _MATELM1( 8 ,8)  += _term;
 /*REACTION*/
  /* ~ i6 <-> is2 ( aS2 , bS )*/
 f_flux =  aS2 * i6 ;
 b_flux =  bS * is2 ;
 _RHS1( 8) -= (f_flux - b_flux);
 _RHS1( 6) += (f_flux - b_flux);
 
 _term =  aS2 ;
 _MATELM1( 8 ,8)  += _term;
 _MATELM1( 6 ,8)  -= _term;
 _term =  bS ;
 _MATELM1( 8 ,6)  -= _term;
 _MATELM1( 6 ,6)  += _term;
 /*REACTION*/
  /* ~ c1 <-> i1 ( Con , Coff )*/
 f_flux =  Con * c1 ;
 b_flux =  Coff * i1 ;
 _RHS1( 5) -= (f_flux - b_flux);
 _RHS1( 13) += (f_flux - b_flux);
 
 _term =  Con ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 13 ,5)  -= _term;
 _term =  Coff ;
 _MATELM1( 5 ,13)  -= _term;
 _MATELM1( 13 ,13)  += _term;
 /*REACTION*/
  /* ~ c2 <-> i2 ( Con * a , Coff / a )*/
 f_flux =  Con * a * c2 ;
 b_flux =  Coff / a * i2 ;
 _RHS1( 4) -= (f_flux - b_flux);
 _RHS1( 12) += (f_flux - b_flux);
 
 _term =  Con * a ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 12 ,4)  -= _term;
 _term =  Coff / a ;
 _MATELM1( 4 ,12)  -= _term;
 _MATELM1( 12 ,12)  += _term;
 /*REACTION*/
  /* ~ c3 <-> i3 ( Con * pow( a , 2.0 ) , Coff / pow( a , 2.0 ) )*/
 f_flux =  Con * pow( a , 2.0 ) * c3 ;
 b_flux =  Coff / pow( a , 2.0 ) * i3 ;
 _RHS1( 3) -= (f_flux - b_flux);
 _RHS1( 11) += (f_flux - b_flux);
 
 _term =  Con * pow( a , 2.0 ) ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 11 ,3)  -= _term;
 _term =  Coff / pow( a , 2.0 ) ;
 _MATELM1( 3 ,11)  -= _term;
 _MATELM1( 11 ,11)  += _term;
 /*REACTION*/
  /* ~ c4 <-> i4 ( Con * pow( a , 3.0 ) , Coff / pow( a , 3.0 ) )*/
 f_flux =  Con * pow( a , 3.0 ) * c4 ;
 b_flux =  Coff / pow( a , 3.0 ) * i4 ;
 _RHS1( 2) -= (f_flux - b_flux);
 _RHS1( 10) += (f_flux - b_flux);
 
 _term =  Con * pow( a , 3.0 ) ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 10 ,2)  -= _term;
 _term =  Coff / pow( a , 3.0 ) ;
 _MATELM1( 2 ,10)  -= _term;
 _MATELM1( 10 ,10)  += _term;
 /*REACTION*/
  /* ~ c5 <-> i5 ( Con * pow( a , 4.0 ) , Coff / pow( a , 4.0 ) )*/
 f_flux =  Con * pow( a , 4.0 ) * c5 ;
 b_flux =  Coff / pow( a , 4.0 ) * i5 ;
 _RHS1( 1) -= (f_flux - b_flux);
 _RHS1( 9) += (f_flux - b_flux);
 
 _term =  Con * pow( a , 4.0 ) ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 9 ,1)  -= _term;
 _term =  Coff / pow( a , 4.0 ) ;
 _MATELM1( 1 ,9)  -= _term;
 _MATELM1( 9 ,9)  += _term;
 /*REACTION*/
  /* ~ o <-> i6 ( Oon , Ooff )*/
 f_flux =  Oon * o ;
 b_flux =  Ooff * i6 ;
 _RHS1( 8) += (f_flux - b_flux);
 
 _term =  Oon ;
 _MATELM1( 8 ,0)  -= _term;
 _term =  Ooff ;
 _MATELM1( 8 ,8)  += _term;
 /*REACTION*/
   /* c1 + c2 + c3 + c4 + c5 + i1 + i2 + i3 + i4 + i5 + i6 + is1 + is2 + o = 1.0 */
 _RHS1(0) =  1.0;
 _MATELM1(0, 0) = 1;
 _RHS1(0) -= o ;
 _MATELM1(0, 6) = 1;
 _RHS1(0) -= is2 ;
 _MATELM1(0, 7) = 1;
 _RHS1(0) -= is1 ;
 _MATELM1(0, 8) = 1;
 _RHS1(0) -= i6 ;
 _MATELM1(0, 9) = 1;
 _RHS1(0) -= i5 ;
 _MATELM1(0, 10) = 1;
 _RHS1(0) -= i4 ;
 _MATELM1(0, 11) = 1;
 _RHS1(0) -= i3 ;
 _MATELM1(0, 12) = 1;
 _RHS1(0) -= i2 ;
 _MATELM1(0, 13) = 1;
 _RHS1(0) -= i1 ;
 _MATELM1(0, 1) = 1;
 _RHS1(0) -= c5 ;
 _MATELM1(0, 2) = 1;
 _RHS1(0) -= c4 ;
 _MATELM1(0, 3) = 1;
 _RHS1(0) -= c3 ;
 _MATELM1(0, 4) = 1;
 _RHS1(0) -= c2 ;
 _MATELM1(0, 5) = 1;
 _RHS1(0) -= c1 ;
 /*CONSERVATION*/
   } return _reset;
 }
 
static int  rates ( _threadargsprotocomma_ double _lv ) {
   alpha = a0 * exp ( ( _lv - vha ) / vca ) ;
   beta = b0 * exp ( ( _lv - vhb ) / vcb ) ;
   gamma = g0 ;
   delta = d0 ;
   a = pow( ( ( Coff / Con ) / ( Ooff / Oon ) ) , ( 1.0 / 8.0 ) ) ;
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
 
double modulationDA ( _threadargsproto_ ) {
   double _lmodulationDA;
 _lmodulationDA = 1.0 + modDA * ( maxModDA - 1.0 ) * levelDA ;
   
return _lmodulationDA;
 }
 
static void _hoc_modulationDA(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  modulationDA ( _p, _ppvar, _thread, _nt );
 hoc_retpushx(_r);
}
 
/*CVODE ode begin*/
 static int _ode_spec1(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
 {int _i; for(_i=0;_i<14;_i++) _p[_dlist1[_i]] = 0.0;}
 rates ( _threadargscomma_ v ) ;
 /* ~ c1 <-> c2 ( 4.0 * alpha , beta )*/
 f_flux =  4.0 * alpha * c1 ;
 b_flux =  beta * c2 ;
 Dc1 -= (f_flux - b_flux);
 Dc2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c2 <-> c3 ( 3.0 * alpha , 2.0 * beta )*/
 f_flux =  3.0 * alpha * c2 ;
 b_flux =  2.0 * beta * c3 ;
 Dc2 -= (f_flux - b_flux);
 Dc3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c3 <-> c4 ( 2.0 * alpha , 3.0 * beta )*/
 f_flux =  2.0 * alpha * c3 ;
 b_flux =  3.0 * beta * c4 ;
 Dc3 -= (f_flux - b_flux);
 Dc4 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c4 <-> c5 ( alpha , 4.0 * beta )*/
 f_flux =  alpha * c4 ;
 b_flux =  4.0 * beta * c5 ;
 Dc4 -= (f_flux - b_flux);
 Dc5 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c5 <-> o ( gamma , delta )*/
 f_flux =  gamma * c5 ;
 b_flux =  delta * o ;
 Dc5 -= (f_flux - b_flux);
 Do += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ o <-> is1 ( aS1 , bS )*/
 f_flux =  aS1 * o ;
 b_flux =  bS * is1 ;
 Do -= (f_flux - b_flux);
 Dis1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ i1 <-> i2 ( 4.0 * alpha * a , beta / a )*/
 f_flux =  4.0 * alpha * a * i1 ;
 b_flux =  beta / a * i2 ;
 Di1 -= (f_flux - b_flux);
 Di2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ i2 <-> i3 ( 3.0 * alpha * a , 2.0 * beta / a )*/
 f_flux =  3.0 * alpha * a * i2 ;
 b_flux =  2.0 * beta / a * i3 ;
 Di2 -= (f_flux - b_flux);
 Di3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ i3 <-> i4 ( 2.0 * alpha * a , 3.0 * beta / a )*/
 f_flux =  2.0 * alpha * a * i3 ;
 b_flux =  3.0 * beta / a * i4 ;
 Di3 -= (f_flux - b_flux);
 Di4 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ i4 <-> i5 ( alpha * a , 4.0 * beta / a )*/
 f_flux =  alpha * a * i4 ;
 b_flux =  4.0 * beta / a * i5 ;
 Di4 -= (f_flux - b_flux);
 Di5 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ i5 <-> i6 ( gamma , delta )*/
 f_flux =  gamma * i5 ;
 b_flux =  delta * i6 ;
 Di5 -= (f_flux - b_flux);
 Di6 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ i6 <-> is2 ( aS2 , bS )*/
 f_flux =  aS2 * i6 ;
 b_flux =  bS * is2 ;
 Di6 -= (f_flux - b_flux);
 Dis2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c1 <-> i1 ( Con , Coff )*/
 f_flux =  Con * c1 ;
 b_flux =  Coff * i1 ;
 Dc1 -= (f_flux - b_flux);
 Di1 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c2 <-> i2 ( Con * a , Coff / a )*/
 f_flux =  Con * a * c2 ;
 b_flux =  Coff / a * i2 ;
 Dc2 -= (f_flux - b_flux);
 Di2 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c3 <-> i3 ( Con * pow( a , 2.0 ) , Coff / pow( a , 2.0 ) )*/
 f_flux =  Con * pow( a , 2.0 ) * c3 ;
 b_flux =  Coff / pow( a , 2.0 ) * i3 ;
 Dc3 -= (f_flux - b_flux);
 Di3 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c4 <-> i4 ( Con * pow( a , 3.0 ) , Coff / pow( a , 3.0 ) )*/
 f_flux =  Con * pow( a , 3.0 ) * c4 ;
 b_flux =  Coff / pow( a , 3.0 ) * i4 ;
 Dc4 -= (f_flux - b_flux);
 Di4 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ c5 <-> i5 ( Con * pow( a , 4.0 ) , Coff / pow( a , 4.0 ) )*/
 f_flux =  Con * pow( a , 4.0 ) * c5 ;
 b_flux =  Coff / pow( a , 4.0 ) * i5 ;
 Dc5 -= (f_flux - b_flux);
 Di5 += (f_flux - b_flux);
 
 /*REACTION*/
  /* ~ o <-> i6 ( Oon , Ooff )*/
 f_flux =  Oon * o ;
 b_flux =  Ooff * i6 ;
 Do -= (f_flux - b_flux);
 Di6 += (f_flux - b_flux);
 
 /*REACTION*/
   /* c1 + c2 + c3 + c4 + c5 + i1 + i2 + i3 + i4 + i5 + i6 + is1 + is2 + o = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE matsol*/
 static int _ode_matsol1(void* _so, double* _rhs, double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset=0;{
 double b_flux, f_flux, _term; int _i;
   b_flux = f_flux = 0.;
 {int _i; double _dt1 = 1.0/dt;
for(_i=0;_i<14;_i++){
  	_RHS1(_i) = _dt1*(_p[_dlist1[_i]]);
	_MATELM1(_i, _i) = _dt1;
      
} }
 rates ( _threadargscomma_ v ) ;
 /* ~ c1 <-> c2 ( 4.0 * alpha , beta )*/
 _term =  4.0 * alpha ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 4 ,5)  -= _term;
 _term =  beta ;
 _MATELM1( 5 ,4)  -= _term;
 _MATELM1( 4 ,4)  += _term;
 /*REACTION*/
  /* ~ c2 <-> c3 ( 3.0 * alpha , 2.0 * beta )*/
 _term =  3.0 * alpha ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 3 ,4)  -= _term;
 _term =  2.0 * beta ;
 _MATELM1( 4 ,3)  -= _term;
 _MATELM1( 3 ,3)  += _term;
 /*REACTION*/
  /* ~ c3 <-> c4 ( 2.0 * alpha , 3.0 * beta )*/
 _term =  2.0 * alpha ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 2 ,3)  -= _term;
 _term =  3.0 * beta ;
 _MATELM1( 3 ,2)  -= _term;
 _MATELM1( 2 ,2)  += _term;
 /*REACTION*/
  /* ~ c4 <-> c5 ( alpha , 4.0 * beta )*/
 _term =  alpha ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 1 ,2)  -= _term;
 _term =  4.0 * beta ;
 _MATELM1( 2 ,1)  -= _term;
 _MATELM1( 1 ,1)  += _term;
 /*REACTION*/
  /* ~ c5 <-> o ( gamma , delta )*/
 _term =  gamma ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 0 ,1)  -= _term;
 _term =  delta ;
 _MATELM1( 1 ,0)  -= _term;
 _MATELM1( 0 ,0)  += _term;
 /*REACTION*/
  /* ~ o <-> is1 ( aS1 , bS )*/
 _term =  aS1 ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 7 ,0)  -= _term;
 _term =  bS ;
 _MATELM1( 0 ,7)  -= _term;
 _MATELM1( 7 ,7)  += _term;
 /*REACTION*/
  /* ~ i1 <-> i2 ( 4.0 * alpha * a , beta / a )*/
 _term =  4.0 * alpha * a ;
 _MATELM1( 13 ,13)  += _term;
 _MATELM1( 12 ,13)  -= _term;
 _term =  beta / a ;
 _MATELM1( 13 ,12)  -= _term;
 _MATELM1( 12 ,12)  += _term;
 /*REACTION*/
  /* ~ i2 <-> i3 ( 3.0 * alpha * a , 2.0 * beta / a )*/
 _term =  3.0 * alpha * a ;
 _MATELM1( 12 ,12)  += _term;
 _MATELM1( 11 ,12)  -= _term;
 _term =  2.0 * beta / a ;
 _MATELM1( 12 ,11)  -= _term;
 _MATELM1( 11 ,11)  += _term;
 /*REACTION*/
  /* ~ i3 <-> i4 ( 2.0 * alpha * a , 3.0 * beta / a )*/
 _term =  2.0 * alpha * a ;
 _MATELM1( 11 ,11)  += _term;
 _MATELM1( 10 ,11)  -= _term;
 _term =  3.0 * beta / a ;
 _MATELM1( 11 ,10)  -= _term;
 _MATELM1( 10 ,10)  += _term;
 /*REACTION*/
  /* ~ i4 <-> i5 ( alpha * a , 4.0 * beta / a )*/
 _term =  alpha * a ;
 _MATELM1( 10 ,10)  += _term;
 _MATELM1( 9 ,10)  -= _term;
 _term =  4.0 * beta / a ;
 _MATELM1( 10 ,9)  -= _term;
 _MATELM1( 9 ,9)  += _term;
 /*REACTION*/
  /* ~ i5 <-> i6 ( gamma , delta )*/
 _term =  gamma ;
 _MATELM1( 9 ,9)  += _term;
 _MATELM1( 8 ,9)  -= _term;
 _term =  delta ;
 _MATELM1( 9 ,8)  -= _term;
 _MATELM1( 8 ,8)  += _term;
 /*REACTION*/
  /* ~ i6 <-> is2 ( aS2 , bS )*/
 _term =  aS2 ;
 _MATELM1( 8 ,8)  += _term;
 _MATELM1( 6 ,8)  -= _term;
 _term =  bS ;
 _MATELM1( 8 ,6)  -= _term;
 _MATELM1( 6 ,6)  += _term;
 /*REACTION*/
  /* ~ c1 <-> i1 ( Con , Coff )*/
 _term =  Con ;
 _MATELM1( 5 ,5)  += _term;
 _MATELM1( 13 ,5)  -= _term;
 _term =  Coff ;
 _MATELM1( 5 ,13)  -= _term;
 _MATELM1( 13 ,13)  += _term;
 /*REACTION*/
  /* ~ c2 <-> i2 ( Con * a , Coff / a )*/
 _term =  Con * a ;
 _MATELM1( 4 ,4)  += _term;
 _MATELM1( 12 ,4)  -= _term;
 _term =  Coff / a ;
 _MATELM1( 4 ,12)  -= _term;
 _MATELM1( 12 ,12)  += _term;
 /*REACTION*/
  /* ~ c3 <-> i3 ( Con * pow( a , 2.0 ) , Coff / pow( a , 2.0 ) )*/
 _term =  Con * pow( a , 2.0 ) ;
 _MATELM1( 3 ,3)  += _term;
 _MATELM1( 11 ,3)  -= _term;
 _term =  Coff / pow( a , 2.0 ) ;
 _MATELM1( 3 ,11)  -= _term;
 _MATELM1( 11 ,11)  += _term;
 /*REACTION*/
  /* ~ c4 <-> i4 ( Con * pow( a , 3.0 ) , Coff / pow( a , 3.0 ) )*/
 _term =  Con * pow( a , 3.0 ) ;
 _MATELM1( 2 ,2)  += _term;
 _MATELM1( 10 ,2)  -= _term;
 _term =  Coff / pow( a , 3.0 ) ;
 _MATELM1( 2 ,10)  -= _term;
 _MATELM1( 10 ,10)  += _term;
 /*REACTION*/
  /* ~ c5 <-> i5 ( Con * pow( a , 4.0 ) , Coff / pow( a , 4.0 ) )*/
 _term =  Con * pow( a , 4.0 ) ;
 _MATELM1( 1 ,1)  += _term;
 _MATELM1( 9 ,1)  -= _term;
 _term =  Coff / pow( a , 4.0 ) ;
 _MATELM1( 1 ,9)  -= _term;
 _MATELM1( 9 ,9)  += _term;
 /*REACTION*/
  /* ~ o <-> i6 ( Oon , Ooff )*/
 _term =  Oon ;
 _MATELM1( 0 ,0)  += _term;
 _MATELM1( 8 ,0)  -= _term;
 _term =  Ooff ;
 _MATELM1( 0 ,8)  -= _term;
 _MATELM1( 8 ,8)  += _term;
 /*REACTION*/
   /* c1 + c2 + c3 + c4 + c5 + i1 + i2 + i3 + i4 + i5 + i6 + is1 + is2 + o = 1.0 */
 /*CONSERVATION*/
   } return _reset;
 }
 
/*CVODE end*/
 
static int _ode_count(int _type){ return 14;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
  }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 14; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _cvode_sparse_thread(&_thread[_cvspth1]._pvoid, 14, _dlist1, _p, _ode_matsol1, _ppvar, _thread, _nt);
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
  ena = _ion_ena;
 _ode_matsol_instance1(_threadargs_);
 }}
 
static void _thread_cleanup(Datum* _thread) {
   _nrn_destroy_sparseobj_thread(_thread[_cvspth1]._pvoid);
   _nrn_destroy_sparseobj_thread(_thread[_spth1]._pvoid);
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  ct = ct0;
  c5 = c50;
  c4 = c40;
  c3 = c30;
  c2 = c20;
  c1 = c10;
  it = it0;
  ist = ist0;
  is2 = is20;
  is1 = is10;
  ift = ift0;
  i6 = i60;
  i5 = i50;
  i4 = i40;
  i3 = i30;
  i2 = i20;
  i1 = i10;
  o = o0;
 {
    _ss_sparse_thread(&_thread[_spth1]._pvoid, 14, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 14; ++_i) {
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
  ena = _ion_ena;
 initmodel(_p, _ppvar, _thread, _nt);
 }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   g = gbar * o * modulationDA ( _threadargs_ ) ;
   ina = g * ( v - ena ) ;
   ct = c1 + c2 + c3 + c4 + c5 ;
   ift = i1 + i2 + i3 + i4 + i5 + i6 ;
   ist = is1 + is2 ;
   it = ift + ist ;
   }
 _current += ina;

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
  ena = _ion_ena;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dina;
  _dina = ina;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
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
  ena = _ion_ena;
 {  sparse_thread(&_thread[_spth1]._pvoid, 14, _slist1, _dlist1, _p, &t, dt, kin, _linmat1, _ppvar, _thread, _nt);
     if (secondorder) {
    int _i;
    for (_i = 0; _i < 14; ++_i) {
      _p[_slist1[_i]] += dt*_p[_dlist1[_i]];
    }}
 } }}
 dt = _dtsav;
}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(o) - _p;  _dlist1[0] = &(Do) - _p;
 _slist1[1] = &(c5) - _p;  _dlist1[1] = &(Dc5) - _p;
 _slist1[2] = &(c4) - _p;  _dlist1[2] = &(Dc4) - _p;
 _slist1[3] = &(c3) - _p;  _dlist1[3] = &(Dc3) - _p;
 _slist1[4] = &(c2) - _p;  _dlist1[4] = &(Dc2) - _p;
 _slist1[5] = &(c1) - _p;  _dlist1[5] = &(Dc1) - _p;
 _slist1[6] = &(is2) - _p;  _dlist1[6] = &(Dis2) - _p;
 _slist1[7] = &(is1) - _p;  _dlist1[7] = &(Dis1) - _p;
 _slist1[8] = &(i6) - _p;  _dlist1[8] = &(Di6) - _p;
 _slist1[9] = &(i5) - _p;  _dlist1[9] = &(Di5) - _p;
 _slist1[10] = &(i4) - _p;  _dlist1[10] = &(Di4) - _p;
 _slist1[11] = &(i3) - _p;  _dlist1[11] = &(Di3) - _p;
 _slist1[12] = &(i2) - _p;  _dlist1[12] = &(Di2) - _p;
 _slist1[13] = &(i1) - _p;  _dlist1[13] = &(Di1) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "na2_ch.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "NA2_CH.MOD\n"
  "\n"
  "c1 - c2 - c3 - c4 - c5 - o  - is1\n"
  "|    |    |    |    |    |\n"
  "i1 - i2 - i3 - i4 - i5 - i6 - is2\n"
  "\n"
  "FAST\n"
  "\n"
  "6/18/2003\n"
  "\n"
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
  "\n"
  "NEURON {\n"
  "	SUFFIX na2_ch\n"
  "	USEION na READ ena WRITE ina\n"
  "	RANGE g, ina, gbar, a\n"
  "	GLOBAL Con, Coff, Oon, Ooff\n"
  "	GLOBAL a0, vha, vca\n"
  "	GLOBAL b0, vhb, vcb\n"
  "	GLOBAL g0\n"
  "	GLOBAL d0\n"
  "	GLOBAL aS1, aS2, bS\n"
  "	RANGE modDA, maxModDA, levelDA\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(mV) = (millivolt)\n"
  "	(mA) = (milliamp)\n"
  "	(S) = (siemens)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	gbar = 1	(S/cm2)\n"
  "\n"
  "	a0 = 37		(1/ms)	: alpha\n"
  "	vha  = 45	(mV)\n"
  "	vca = 40	(mV)\n"
  "\n"
  "	b0 = 10		(1/ms)	: beta\n"
  "	vhb = -50	(mV)\n"
  "	vcb = -10	(mV)\n"
  "\n"
  "	g0 = 40		(1/ms)	: gamma\n"
  "\n"
  "	d0 = 30		(1/ms)	: delta\n"
  "\n"
  "	aS1 = 0.0025	(1/ms)\n"
  "	aS2 = 0.0002	(1/ms)\n"
  "	bS = 0.00017	(1/ms)\n"
  "\n"
  "	Con = 0.001	(1/ms)\n"
  "	Coff = 0.1	(1/ms)\n"
  "	Oon = 1.6	(1/ms)\n"
  "	Ooff = 0.01	(1/ms)\n"
  "        modDA = 0\n"
  "        maxModDA = 1\n"
  "        levelDA = 0\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v	(mV)\n"
  "	ena	(mV)\n"
  "	g	(S/cm2)\n"
  "	ina	(mA/cm2)\n"
  "	alpha	(1/ms)\n"
  "	beta	(1/ms)\n"
  "	gamma	(1/ms)\n"
  "	delta	(1/ms)\n"
  "	a\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	c1  : closed\n"
  "	c2\n"
  "	c3\n"
  "	c4\n"
  "	c5\n"
  "	ct  : total closed\n"
  "	o   : open\n"
  "	i1  : fast inactivated\n"
  "	i2\n"
  "	i3\n"
  "	i4\n"
  "	i5\n"
  "	i6   \n"
  "	ift : total fast inactivated\n"
  "	is1 : slow inactivated\n"
  "	is2\n"
  "	ist : total slow inactivated\n"
  "	it  : total inactivated\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE kin METHOD sparse\n"
  "	g = gbar*o*modulationDA()\n"
  "	ina = g*(v-ena)\n"
  "	ct = c1 + c2 + c3 + c4 + c5\n"
  "	ift = i1 + i2 + i3 + i4 + i5 + i6\n"
  "	ist = is1 + is2\n"
  "	it = ift + ist\n"
  "}\n"
  "\n"
  "INITIAL {\n"
  "	SOLVE kin STEADYSTATE sparse\n"
  "}\n"
  "\n"
  "KINETIC kin{\n"
  "	rates(v)\n"
  "\n"
  "	~ c1 <-> c2 (4*alpha, beta)\n"
  "	~ c2 <-> c3 (3*alpha, 2*beta)\n"
  "	~ c3 <-> c4 (2*alpha, 3*beta)\n"
  "	~ c4 <-> c5 (alpha, 4*beta)\n"
  "	~ c5 <-> o  (gamma, delta)\n"
  "	~ o <-> is1 (aS1, bS)\n"
  "\n"
  "	~ i1 <-> i2 (4*alpha*a, beta/a)\n"
  "	~ i2 <-> i3 (3*alpha*a, 2*beta/a)\n"
  "	~ i3 <-> i4 (2*alpha*a, 3*beta/a)\n"
  "	~ i4 <-> i5 (alpha*a, 4*beta/a)\n"
  "	~ i5 <-> i6 (gamma, delta)\n"
  "	~ i6 <-> is2 (aS2, bS)\n"
  "\n"
  "	~ c1 <-> i1 (Con, Coff)\n"
  "	~ c2 <-> i2 (Con*a, Coff/a)\n"
  "	~ c3 <-> i3 (Con*a^2, Coff/a^2)\n"
  "	~ c4 <-> i4 (Con*a^3, Coff/a^3)\n"
  "	~ c5 <-> i5 (Con*a^4, Coff/a^4)\n"
  "	~ o <-> i6  (Oon, Ooff)\n"
  "\n"
  "	CONSERVE c1+c2+c3+c4+c5+i1+i2+i3+i4+i5+i6+is1+is2+o=1\n"
  "}\n"
  "\n"
  "PROCEDURE rates(v(millivolt)) {\n"
  "	alpha = a0*exp((v-vha)/vca)\n"
  "	beta = b0*exp((v-vhb)/vcb)\n"
  "	gamma = g0\n"
  "	delta = d0\n"
  "\n"
  "	a = ((Coff/Con)/(Ooff/Oon))^(1/8)\n"
  "}\n"
  "\n"
  "FUNCTION modulationDA() {\n"
  "    : returns modulation factor\n"
  "    \n"
  "    modulationDA = 1 + modDA*(maxModDA-1)*levelDA \n"
  "}\n"
  ;
#endif
