$Title  MPB spread model with the treatments of the infested sites - July 2023

$Ontext
Pest invasion spread with the treatments of newly infested sites -
adapted to model a coarse-scale MPB long-distance spread .
Originally created by D.Yemshanov, August 2023, Adapted by E. Hudgins

$Offtext

$onecho > %outpath%/gurobi.prm
LogFile %outpath%/gurobi1.log
Threads 96 * adjust to number of cores
GURO_PAR_PROCGROUPS 2
nodefilestart 1
softmemlimit 800
$offecho

$onecho > gurobi.opt
TimeLimit  125000
readparams %outpath%/gurobi.prm
MIPGap  0.0005
mipstart 1
mipfocus 2
varhint 1
presolve 2
RINS 200
*method 2
nodefilestart 1
softmemlimit 800
PreSparsify 2
quad 1
iis 1

Threads  16
$offecho

$eolcom //

$offlisting

* This stops the echo print of the input file
Option limcol = 0;
* This stops the print of the column listing
Option limrow = 0;
* This stops the print of the equation listing
$offsymxref
$offsymlist
*This stops the print of a complete cross-reference list of symbols
Option solprint = off;
Option sysout = on;
Option profile = 1;

* -- setting the maximum spread time horizon
$if not set t_max $set t_max 12
* -- setting an initialization time horizon min t_set = 2, max t_set = t_ periods
$if not set t_set $set t_set 2
* -- number ot time steps with the treatmemts but solving the spread for t_ steps
$if not set tstep $set tstep 2
* -- number ot time steps to output in results1.txt and popul_nt.txt
$if not set t_out $set t_out 11

* Input folder/file prefix
$if not set outpath $set outpath "Results_tst"
* Output folder/file prefix
$if not set inpath $set inpath "Input_tst/corr"

* Output header with the parameter names in results1.txt (1=yes, 0=no)
$if not set Pnames $set Pnames 1
* fixed cost portion
$if not set costfix $set costfix 0.7
* variable cost portion per population unit
$if not set costpop $set costpop 20
* treatment efficiency
$if not set eff $set eff 0.63
* population growth rate
$if not set Grate $set Grate 1.6
* treatment budget
$if not set budg $set budg 15
* maximum spread distance, meters
$if not set Mdist $set Mdist 70000
* dispersal kernel adjustment coefficient
$if not set Kcoef $set Kcoef 0.0002
* minimum population density for recently infested sites
$if not set popmin $set popmin 0.004
* population density threshold when the infestation becomes detectable and treatable
$if not set popdet $set popdet 0.008
* population density threshold when the infestation starts spreading propagules to other sites
$if not set popspr $set popspr 0.02
* max population density threshold when the spread stops (population collapse)
$if not set popmax $set popmax 1
* selecting the model type: --mNo=2 - full model, =3-spread only, =1 initialization =5- initialization min cost
$if not set mNo $set mNo 3
* scaling coefficient for the undetected infestation component in the objective equation
$if not set NDspr $set NDspr 0
* scaling coefficient for minimizing the infested area including the detecteable infestation
$if not set Dspr $set Dspr 1

* generating the spread pattern from equation (usually with the spread-only model mNo=3
* once the spread pattern is generated and stored in %inpath% folder the treatment solves
* must load the spread pattern from file %inpath%spread.txt
* m_no=3 - spread only model, no treatments
* m_no=4 - initialization p_inf2 model - min detected inf. area * priority for one period, spread and treatments for one period

* m_no=1 - initialization p_inf1 model - min (detected inf. area * priority) for one period, spread for T periods
* m_no=2 - spread and treatments model - min (detected inf. area * priority). s.t. budget constraints  over T periods

* m_no=5 - initialization p_inf1_cost model - min (cost) + detectd spread area penalty s.t. StS target
* m_no=6 - spread and treatments model - min (cost) + detected spread area penalty s.t. StS target

$if not set genspr $set genspr 0
* file with the spread pattern
$if not set sprfile $set sprfile "spread.txt"
* slow-the-spread target, proportion of new detectable infestations in the treatments vs spread-only scenario
* --StS=0.6 assumes that the treatments reduce the infested detctable area * priorty to 60% from the area in the spread only scenario
$if not set StS $set StS 0.6

* buffer distance (meters) around the infested sites in the spread-only scenario - to generate the spread mask
$if not set buff_msk $set buff_msk 25000


Set t_   time periods - full planning horizon
/
$include %inpath%/t%t_max%.txt
/
Set t(t_) time periods - initialization timespan - the full model is run with t=t_
/
$include %inpath%/t%t_set%.txt
/
Set n   sites
/
$include %inpath%/nodes.txt
/

;
Alias (n,m,k)
;

Parameter
t_id(t_) time period integer identifier
/
$ondelim
$include %inpath%/t%t_max%_ID.txt
$offdelim
/

infest(n)  sites that are infested in period 1

x_coord(n)    x coordinate of site n
/
$ondelim
$include %inpath%/x_coord.txt
$offdelim
/
y_coord(n)    y coordinate of site n
/
$ondelim
$include %inpath%/y_coord.txt
$offdelim
/
w_inf0(n)   population density in sites infested in period 1
/
$ondelim
$include %inpath%/infest_w.txt
$offdelim
/
prior(n)    relative priority of slowing the spread in site n
/
$ondelim
$include %inpath%/prior.txt
$offdelim
/
host(n)    host density in site n [0-1]
/
$ondelim
$include %inpath%/host.txt
$offdelim
/
G(n)    growth rate in site n [1-5]
/
$ondelim
$include %inpath%/g_2024.txt
$offdelim
/
w_spr(n)    spread phase threshold
/
$ondelim
$include %inpath%/w_spr_2024.txt
$offdelim
/

* -- comment $onwarning and $offwarning to enable the check for domain violation errors
$onwarning

max_obj_t(t)    objective value (number of inf sites * priority) in period t in the spread-only scenario
/
$ondelim
$include %outpath%/max_obj_t.txt
$offdelim
/

spread(n,n,t_) spread between pairs of sites n in period t
/
$ondelim
$include %inpath%/%sprfile%
$offdelim
/

mask(n,t_)     mask of eligible sites n based on the spread-only scenario.
/
$ondelim
$include %outpath%/mask_spread.txt
$offdelim
/

* --- loading decision variables from $include file
*$include load_other_variables.gms
$include load_all_variables.gms

* alternatively, comment the include file above and uncomment the text below to avoid initialization

$ontext
* - Initializing other decision variables without loading from the files -------
u0(n,t)   u - site n is infested in period t
q0(n,t)   q - treatments of site n in period t
w0(n,t)     w - population density in site n in period t
z0(n,t)     z - site n can spread propagules in period t
v0(n,t)     v - infestation in site n is detectable and treatable in period t
x0(n,t)     x - first infestation of site in period t over timespan t=2...T
d0(n,t)     d - population density in site n in period t is >= w_spr
a0(n,t)     a - population density in site n in period t is <= w_max
zx0(n,t)    zx -  product of x(n_t) and z(n_t)
xq10(n,t)   xq1 - product of x(n_t) and q(n_t)
w10(n,t)    w1 -  expected pop density in site n in period t just before treatments
wx0(n,t)    wx -  product of w(n_t-1) and x(n_t)
wq0(n,t)    wq -  product of w(n_t) and 1-q(n_t)
wq10(n,t)   wq1 - product of w(n_t-1) and q(n_t)
wq20(n,t)   wq2 - product of w(n_t) and q(n_t)
wxq10(n,t)  wxq1- product of wx(n_t) and q(n_t)
$offtext
;

$offwarning
;



Scalar
m_no       model ID - taken from command line parameter --mNo
NZ         number of non-zeroes in the model
P_names    switch to output the header with the parameter names in results1.txt
G0         population growth rate in period t vs t-1
B          upper bound on the treatment budget in period t
Max_dist   maximum spread distance from site n (limits the spread kernel radius)
Kern_coef  dispersal kernel scaling factor (for calibration)
e          treatment efficiency e=]0-1[ - better to stay within 0.15-0.85 to avoid infeasibilities
cost_fix   fixed treatment cost portion for site n in period t
cost_pop   variable treament cost portion - linearly dependent on the population density in site n in period t

w_min      minimum population density in site n at the time of infestation (assigned when x(n_t) = 1)
w_det      minimum population density in site n that enables detection and treatments
w_spr0     baseline minimum population density in site n that enables spread
w_max      maximum population density in site n that stops spread and triggers the population collapse
t_step     number of periods to predict the treatments while predicting the spread for t_set time periods
Sts_target proportional target reduction of the infested area in period t vs the spread-only scenario
gen_spread option to generate the spread pattern and store it in %inpath%spread.txt file
buffer_mask buffer distance (meters) around the infested sites in the spread-only scenario - to generate the spread mask

;
buffer_mask = %buff_msk%; //25000;  // buffer distances from the infestd sites to calculate the spread mask
m_no = %mNo%; // model number - taken from command line parameters
gen_spread = %genspr%;
StS_target = %StS%;   // target proportion of detected infested cells vs spread-only scenario
t_step = %tstep%;
P_names = %Pnames%;
cost_fix = %costfix%; // fixed cost portion to treat site n in period t
cost_pop = %costpop%; // variable cost portion to treat one population unit in site n in period t
e = %eff%;            // treatment+detection efficiency ]0;1[
G0 = %Grate%;         // population growth rate in period t vs. period t-1
B = %budg%;           // an upper bound on the treatment budget
Max_dist = %Mdist%;   // maximum spread distance (meters) from site n (limits the dispersal distance from n)
Kern_coef = %Kcoef%;  // 0.0001; // calibration coefficient in the spread kernel equation

w_min = %popmin%;     // minimum population density that is assigned at the time of infestation
w_det = %popdet%;     // minimum population density that enables detection and treatments
w_spr0 = %popspr%;    // minimum population density that enables spread to other sites
w_max = %popmax%;     // maximum population density at which the population can spread propagules to other sites

* --- reset mask(n,t) to 1 for all sites in the spread-only scenario (m_No = 3)
mask(n,t)$(m_No eq 3) = 1;

*Parameter w_spr(n), G(n);   // Generating the site-specific spread rate and min spread density based on host abundance
//w_spr(n)= w_spr0; // start with the baseline min spread density w_spr0
////w_spr(n)$(host(n) < 0.12 and host(n) ge 0.06) = w_spr0*1.5;
//w_spr(n)$(host(n) < 0.11 ) = w_spr0*2; // increase the min spread density threshold for sites with low host abundance

*w_spr(n)$(w_spr(n) > 0.026) = 0.026;
*w_spr(n)$(w_spr(n) < 0.014) = 0.014;
*//w_spr(n) = w_spr0;
*//G(n)=g0;
*G(n)$(G(n) >1.885)=1.885;
*G(n)$(G(n) <1.015)=1.015;


w_inf0(n)$(w_inf0(n)>w_min * 2) = w_inf0(n) / 2; // conditioning the high infestation densities for previously infested sites at t=1
infest(n)$(w_inf0(n)>0) = 1;  // all sites with w_inf(0)>0 in period 1 are considered infested

* calculating the randomized spread pattern for sites n ------------------------
* spread(n,n,t) is a binary spread indicator between a pair of sites n in period t,
Parameter coin(n,n,t), P_value(n,n,t), K_coef(n,n,t), dist(n,n), spread0(n);

* Kernel equation for MPB (provided by Emma Hudgins, July 2023):

dist(n,m) = sqrt(
(x_coord(n) - x_coord(m)) * (x_coord(n) - x_coord(m)) +
(y_coord(n) - y_coord(m)) * (y_coord(n) - y_coord(m)) ); // distance between sites n and m


* -- generating the spread patterns - spread(n,m,t)
* -- use swith: --genspr=1 to generate the spread pattern and store it in spread.txt,
* -- otherwise: load the spread pattern from %inpath%/spread.txt
file fw1 / '%inpath%/spread.txt' /;   // file with the spread patterns
fw1.pc = 6; fw1.nd = 0; fw1.pw=32767;
scalar kernel_coef;   kernel_coef = -0.5615539;
if(gen_spread eq 1,
   coin(n,m,t) = uniform(0,1); // draw a uniform random number from interval [0;1]
   K_coef(n,m,t)$(dist(n,m) le Max_dist) = exp(kernel_coef * dist(n,m) * Kern_coef); //P(spead) =f(dist)
   spread(n,m,t) = 0;                                // reset the spread pattern
   spread(n,m,t)$(coin(n,m,t) le K_coef(n,m,t)) = 1; // generate the binary spread patterns
   loop(n,loop(m,loop(t$(spread(n,m,t)>0),put fw1 n.tl m.tl t.tl spread(n,m,t) /)));
   );
putclose fw1;

Scalar M_0, M_; // big-M value for the population density constraints
M_0 = smax(n, w_inf0(n)) * G0;
loop(t_, M_ = M_0 * G0; M_0 = M_); M_ = M_ * 5;

* ------------------------------------------------------------------------------
Variable
   obj_val    objective value
Binary variable
   q(n,t)     site n is treated in period t
   u(n,t)     site n is infested in period t
   v(n,t)     site n is detectable and can be treated in period t
   d(n,t)     population density in site n in period t is >= w_spr
   a(n,t)     population density in site n in period t is <= w_max
   z(n,t)     site n can spread propagules in period t elsewhere
   x(n,t)     first infestation indicator in site n in period t
   zx(n,t)    product of x(n_t) and z(n_t)
   xq(n,t)    product of x(n_t) and 1-q(n_t)
   xq1(n,t)   product of x(n_t) and q(n_t)

Positive variable
   w(n,t)     population density in site n in period t
   w1(n,t)    expected population density in site n in period t just before the treatments
   wx(n,t)    product of w(n_t-1) and x(n_t)
   wq(n,t)    product of w(n_t) and 1-q(n_t)
   wq1(n,t)   product of w(n_t-1) and q(n_t)
   wq2(n,t)   product of w(n_t) and q(n_t)
   w1q(n,t)   product of w1(n_t) and 1-q(n_t)
   wxq(n,t)   product of wx(n_t) and 1-q(n_t)
   wxq1(n,t)  product of wx(n_t) and q(n_t)
   V_budg(t)  penalty for unspent budget below B * 0.95
   B_(t)      treatment budget for period t
;

* -- Initializing the decision variables ---------------------------------------
spread(n,m,t)$(spread(n,m,t) < 1 and mask(n,t)<1) = no;
w.fx(n,"1") = w_inf0(n);  // initializing the population density in period 1
u.fx(n,"1") = infest(n);  // initializing the infested sites in period 1
w1.fx(n,"1")$(w_inf0(n)>0) = w_inf0(n); // initializing the population density in period 1
x.fx(n,"1")$(infest(n)>0) = 1;


* -- initializing the decision variables from the loaded values ----------------
u.l(n,t) = u0(n,t);       w.l(n,t)  = w0(n,t);
x.l(n,t) = x0(n,t);       w1.l(n,t) = w10(n,t);
v.l(n,t) = v0(n,t);       z.l(n,t)  = z0(n,t);      xq1.l(n,t) = xq10(n,t);
d.l(n,t) = d0(n,t);       a.l(n,t)  = a0(n,t);      q.l(n,t)   = q0(n,t);
wxq1.l(n,t) = wxq10(n,t); wx.l(n,t) = wx0(n,t);     zx.l(n,t)  = zx0(n,t);
wq1.l(n,t) = wq10(n,t);   wq.l(n,t) = wq0(n,t);     wq2.l(n,t) = wq20(n,t);

$ontext
u.l(n,t) = 0;       w.l(n,t) = 0;       q.l(n,t) = 0;      v.l(n,t) = 0;
z.l(n,t) = 0;       x.l(n,t) = 0;       d.l(n,t) = 0;      a.l(n,t) = 0;
xq1.l(n,t) = 0;     w1.l(n,t) = 0;      wxq1.l(n,t) = 0;   wx.l(n,t) = 0;
zx.l(n,t) = 0;      wq1.l(n,t) = 0;     wq.l(n,t) = 0;     wq2.l(n,t) = 0;
$offtext

obj_val.l = 0;

u.fx(n,t)$(mask(n,t)<1) = 0; v.fx(n,t)$(mask(n,t)<1) = 0;

if (m_no eq 1 or m_no eq 2 or m_no eq 4, q.fx(n,t)$(t_ID(t) > t_step) =  0;);
if (m_no eq 3, q.fx(n,t) =0;);

* ------------------------------------------------------------------------------

Equations
  obj_N             objective value - min N(infested sites) times priority over t periods
  obj_N_x           objective value - min N(first infestations) times priority over t periods
  obj_N_v           objective value - min N(detectable infestations) times priority over t periods
  obj_N_v1          objective value - min N(detectable infestations) times priority over t_step periods
  obj_N_u           objective value - min N(infsted sites) at time t = T (card(t))
  obj_cost          minimizing the treatment cost
  budget_n(t)       max treatment budget in period t - fixed site treatment cost = 1
  budget(t)         max treatment budget in period t - population-dependent cost of treating site n
  budgetG(t)        max treatment budget in period t - population-dependent cost of treating site n
  V_budgetG(t)      penalty for the unspent budget below 0.95 * B
  budget_min(t)     min treatment budget in period t - population-dependent cost of treating site n
  buff_nodes(n,n,t) spread to node n from other infested nodes in period t
  buff_wmin(n,n,t)  newly infested nodes are assigned the population density w_min
  sum_buff_n(n,t)   x(n_t) = 0 if there is no new spread to n in period t from other sites
  u_inf_nt(n,t)     site n is infested in period t if w_nt >= w_min
  u_uninf_nt(n,t)   site n is not infested in period t if w_nt <= w(nt) div w_min
  w_inf_u(n,t)      a newly infested site needs to have population density >= wmin

  pop_nt(n,t)       population growth in site n in period t vs period t-1
  pop_nt1(n,t)      population growth in site n in period t vs period t-1 (v2)
  pop_nt1_(n,t)     population growth in site n in period t vs period t-1 (v3 simplified)
  pop_nt0(n,t)      population growth in site n in period t vs period t-1 expected w-o treatments
  pop_nt00(n,t)     population growth in site n in period t vs period t-1 spread only
  pop_nt00_(n,t)    population growth in site n in period t vs period t-1 spread only

  w_det1(n,t)       site is detectable-treatable if w >= w_det
  w_det2(n,t)       site is detectable-treatable if w >= w_det
  w_wspr(n,t)       w(n_t) - w_spr =l= M_ * a(n_t)
  w_wspr1(n,t)      w(n_t) =g= w_spr * z(n_t)
  w_wspr2(n,t)      z(n_t) * w_spr =l= w(n_t)
  w_wmax(n,t)       w_max - w(n_t) =l= M_ * d(n_t)
  w_wmax1(n,t)      w(n_t) =l= w_max * z(n_t) + M_* (1-z(n_t))

  z_a_d1(n,t)       z(n_t) =l= a(n_t)
  z_a_d2(n,t)       z(n_t) =l= d(n_t)
  z_a_d3(n,t)       z(n_t) =g= a(n_t) + d(n_t) - 1
  x_u(n,t)          x(n_t) =g= u(n_t) - u(n_t-1)
  x_u1(n,t)         x(n_t) + u(n_t) + u(n_t-1) =l= 2
  x_u2(n,t)         x(n_t) =l= u(n_t)
  x_u3(n,t)         x(n_t) =l= 1-u(n_t-1)
  wx_1(n,t)         wx(n_t) =l= M_ * x(n_t)
  wx_2(n,t)         w(n_t-1) - wx(n_t) - M_ * (1-x(n_t)) =l= 0
  wx_3(n,t)         w(n_t) =g= wx(n_t)
  zx_1(n,t)         product of z(n_t) and x(n_t)
  zx_2(n,t)         product of z(n_t) and x(n_t)
  zx_3(n,t)         product of z(n_t) and x(n_t)
  wxq_1(n,t)        product of wx(n_t) and 1-q(n_t)
  wxq_2(n,t)        product of wx(n_t) and 1-q(n_t)
  wxq_3(n,t)        product of wx(n_t) and 1-q(n_t)
  wxq1_1(n,t)       product of wx(n_t) and q(n_t)
  wxq1_2(n,t)       product of wx(n_t) and q(n_t)
  wxq1_3(n,t)       product of wx(n_t) and q(n_t)
  wq_1(n,t)         product of w(n_t-1) and 1-q(n_t)
  wq_2(n,t)         product of w(n_t-1) and 1-q(n_t)
  wq_3(n,t)         product of w(n_t-1) and 1-q(n_t)
  w1q_1(n,t)        product of w1(n_t) and 1-q(n_t)
  w1q_2(n,t)        product of w1(n_t) and 1-q(n_t)
  w1q_3(n,t)        product of w1(n_t) and 1-q(n_t)
  wq1_1(n,t)        product of w(n_t-1) and q(n_t)
  wq1_2(n,t)        product of w(n_t-1) and q(n_t)
  wq1_3(n,t)        product of w(n_t-1) and q(n_t)
  wq2_1(n,t)        product of w(n_t) and q(n_t)
  wq2_2(n,t)        product of w(n_t) and q(n_t)
  wq2_3(n,t)        product of w(n_t) and q(n_t)

  u_t(n,t)          w(n_t) =g= wq(n_t)
  u_t0(n,t)         w(n_t) =g= w(n_t-1)
  q_v(n,t)          no treatment of sites with w(n_t) < w_det
  xq_1(n,t)         product of x(n_t) and 1-q(n_t)
  xq_2(n,t)         product of x(n_t) and 1-q(n_t)
  xq_3(n,t)         product of x(n_t) and 1-q(n_t)
  xq1_1(n,t)        product of x(n_t) and q(n_t)
  xq1_2(n,t)        product of x(n_t) and q(n_t)
  xq1_3(n,t)        product of x(n_t) and q(n_t)

  max_N_inf         limit the number of infested sites over t_step periods
  min_N_inf(t)      limit the number of infested sites in period t
  sum_q(n)          site n can only be treated once over T periods
  q_seq(n,t)        no treatments in site n in two consecutive periods t-1 and t

  fix_u_t(n,t)      fix the infested status for periods t-1 when running initialization
  fix_q_t(n,t)      fix the treated status for periods t-1 when running initialization
  fix_w_t(n,t)      fix the pop.densities for periods t-1 when running initialization
  fix_u_t_(n,t)     fix the infested status for periods T-1 when running initialization
  fix_q_t_(n,t)     fix the treated status for periods T-1 when running initialization
  fix_w_t_(n,t)     fix the pop.densities for periods T-1 when running initialization
  fix_q_t_0(n,t)    no treatments in periods > tstep
  fix_q(n,t)        fix the treatment variable q(n_t)
  no_treat(n,t)     q(n_t) =e= 0 - no treatments - spread only
  inf_u_t(n,t)      site n once infested stays infested

;
B_.l(t) = 5e3; // * t_step; // set the initial B_ value sufficiently high
Parameter prior_t(t);   prior_t(t) = 1;// temporal priority in the objective function
if (m_No eq 4, prior_t(t)$(t_ID(t) le t_step) = 1; prior_t(t)$(t_ID(t) > t_step) = 0;);

Scalar ND_spr, D_spr; ND_spr = %NDspr%; D_spr=%Dspr%;
* -- Equations -----------------------------------------------------------------
obj_N..    sum(t$(t_ID(t)>1), sum(n$(mask(n,t)>0), prior(n) * (u(n,t) - u(n,t-1)) ))  =e= obj_val; //min total N infested sites time priority
obj_N_u..  sum(t$(t_ID(t) eq card(t)), sum(n$(mask(n,t)>0), prior(n) * u(n,t) ))  =e= obj_val; //min total N infested sites time priority

* -- main objective ------------------------------------------------------------
obj_N_v..  D_spr * sum(n$(mask(n,"2")>0), prior(n) * v(n,"2"))
         + D_spr * sum(t$(t_ID(t)>2), sum(n$(mask(n,t)>0), prior_t(t) * prior(n) * (v(n,t) - v(n,t-1)) ))
         + ND_spr * sum(t$(t_ID(t)>1), sum(n$(mask(n,t)>0), prior_t(t) * prior(n) * (x(n,t) ) ))
        =e= obj_val; //min N first infestations times priority

* -- min detected infested area * priority only for t <= t_step periods --------
obj_N_v1.. D_spr * sum(n$(mask(n,"2")>0), prior(n) * v(n,"2"))
         + D_spr * sum(t$(t_ID(t)>2 and t_ID(t) le t_step), sum(n$(mask(n,t)>0), prior(n) * (v(n,t) - v(n,t-1)) ))
         + ND_spr * sum(t$(t_ID(t)>1 and t_ID(t) le t_step), sum(n$(mask(n,t)>0), prior(n) * (u(n,t) - u(n,t-1)) ))
        =e= obj_val; //min N first infestations times priority

obj_N_x..  sum(t, sum(n$(mask(n,t)>0), prior(n) * x(n,t) ))  =e= obj_val; //min N first infestations times priority

* -- minimiing the treatment cost over t periods -------------------------------
obj_cost..    D_spr * sum(n$(mask(n,"2")>0), prior(n) * v(n,"2"))
            + D_spr * sum(t$(t_ID(t)>2), sum(n$(mask(n,t)>0), prior(n) * (v(n,t) - v(n,t-1)) ))
            + ND_spr * sum(t$(t_ID(t)>1), sum(n$(mask(n,t)>0), prior(n) * (x(n,t) ) ))
             + 30 * sum(t$(t_ID(t)>1), B_(t))
        =e= obj_val;


* -- treatment budget constraints -----
budget_n(t)$(t_ID(t)>1).. sum(n$(mask(n,t)>0), q(n,t)) =l= B;     // fixed cost=1 of treating site n in period t
budgetG(t)$(t_ID(t)>1)..  sum(n$(mask(n,t)>0), cost_fix * q(n,t) + cost_pop * wq1(n,t) * G(n)) =l= B;   // maximum treatment budget - the variable cost portion depends on the pop. density in site n in period t just before treatments

V_budgetG(t)$(t_ID(t)>1 and t_ID(t) le t_step)..
       sum(n$(mask(n,t)>0), cost_fix * q(n,t) + cost_pop * wq1(n,t) * G(n)) =l= B_(t);   // maximum treatment budget - the variable cost portion depends on the pop. density in site n in period t just before treatments

budget_min(t)$(t_ID(t)>1 and t_ID(t) le t_step)..
       sum(n$(mask(n,t)>0), cost_fix * q(n,t) + cost_pop * wq1(n,t) * G(n)) =g= B_(t) * 0.8;

* -- limit the number of infested sites over t_step periods --------------------
max_N_inf..     sum(n$(mask(n,"2")>0), prior(n) * v(n,"2")) +
                sum(t$(t_ID(t) le t_step and t_ID(t)>2),
                sum(n$(mask(n,t)>0), prior(n) *  (v(n,t) - v(n,t-1))   ))
            =l= StS_target *
                sum(t$(t_ID(t) le t_step and t_ID(t)>1),max_obj_t(t) );
* -- minimum infestation StS target in period t - used for initialization only
min_N_inf(t)$(t_ID(t)>1 and t_ID(t) le t_step)..
//max_N_inf..  sum(t$(t_ID(t)>1 and t_ID(t) le t_step),
             sum(n$(mask(n,t)>0), prior(n) * (v(n,t)-v(n,t-1) )) =g= StS_target * 0.9 *
             max_obj_t(t);

* -- detection / treatment threshold ------
w_det1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   w(n,t) - w_det  =l= v(n,t) * M_ * 5;
w_det2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   v(n,t) =l= w(n,t) / w_det;

* -- z(n,t) = 1 if w_spr <= w(n,t) <= w_max : propagule spread indicator z(n,t) for site n period t
w_wmax(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   w_max - w(n,t) =l= 5* M_ * d(n,t);
w_wmax1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =l= w_max * z(n,t) + M_* (1-z(n,t));
w_wspr(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   w(n,t) - w_spr(n) =l= M_ * a(n,t);
w_wspr1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =g= w_spr(n) * z(n,t);
w_wspr2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  z(n,t) * w_spr(n) =l= w(n,t);
z_a_d1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   z(n,t) =l= a(n,t);
z_a_d2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   z(n,t) =l= d(n,t);
z_a_d3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   z(n,t) =g= a(n,t) + d(n,t) -1 ;


* -- u(n,t)=1 when w(n,t)>0, otherwise =0 - the infested site indicator u(n,t)
u_inf_nt(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    u(n,t) * M_ =g= w(n,t); //  site n is infested in period t if w(n,t) >= w_min
u_uninf_nt(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  u(n,t) =l= w(n,t) * 8000; // site n is not infested in period t if w(n,t) = 0

* -- pop. density in period t >= pop.dens. in period t-1 for the untreated sites n
u_t(n,t)$(t_ID(t)>1 and mask(n,t)>0)..      w(n,t) =g= wq(n,t);
* -- pop. density in period t >= pop.dens. in period t-1 - spread-only scenario
u_t0(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     w(n,t) =g= w(n,t-1);
* -- pop. density in period t >= w_min in period t if site n was infested in period t
w_inf_u(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =g= x(n,t) * w_min;
* -- site n can only be treated once over T peiods ---
sum_q(n)..                  sum(t$(mask(n,t)>0), q(n,t)) =l= 1;
* -- no treatments in site n in two consecutive periods t and t-1
q_seq(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    q(n,t) + q(n,t-1) =l= 1;

* -- product of w(n,t-1) and 1-q(n,t) - pop.density in peiord t-1 and NO treatments in period t
wq_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     wq(n,t) =l= M_ * (1-q(n,t));
wq_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     w(n,t-1) - wq(n,t) - M_ * (q(n,t)) =l= 0;
wq_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     wq(n,t) =l= w(n,t-1);
* -- product of w(n,t-1) and q(n,t) - pop.density in period t-1 and site n is treated in period t
wq1_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wq1(n,t) =l= M_ * (q(n,t));
wq1_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    w(n,t-1) - wq1(n,t) - M_ * (1-q(n,t)) =l= 0;
wq1_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wq1(n,t) =l= w(n,t-1);
* -- product of w(n,t) and q(n,t) - pop.density in period t and site n is treated in period t
wq2_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wq2(n,t) =l= M_ * (q(n,t));
wq2_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    w(n,t) - wq2(n,t) - M_ * (1-q(n,t)) =l= 0;
wq2_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wq2(n,t) =l= w(n,t);
* -- product of x(n,t) and w(n,t-1) --------------------------
wx_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     wx(n,t) =l= M_ * x(n,t);
wx_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     w(n,t-1) - wx(n,t) - M_ * (1-x(n,t)) =l= 0;
wx_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     w(n,t-1) =g= wx(n,t);
* -- product of wx(n,t) and 1-q(n,t) -------------------------
wxq_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wxq(n,t) =l= M_ * (1-q(n,t));
wxq_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wx(n,t) - wxq(n,t) - M_ * (q(n,t)) =l= 0;
wxq_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    wxq(n,t) =l= wx(n,t);
* -- product of wx(n,t) and q(n,t) ---------------------------
wxq1_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   wxq1(n,t) =l= M_ * (q(n,t));
wxq1_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   wx(n,t) - wxq1(n,t) - M_ * (1-q(n,t)) =l= 0;
wxq1_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   wxq1(n,t) =l= wx(n,t);
* -- product of x(n,t) and 1-q(n,t) --------------------------
xq_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     xq(n,t) =l= x(n,t);
xq_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     xq(n,t) =l= (1-q(n,t));
xq_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     xq(n,t) =g= x(n,t) + (1-q(n,t)) - 1;
* -- product of x(n,t) and q(n,t) ----------------------------
xq1_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    xq1(n,t) =l= x(n,t);
xq1_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    xq1(n,t) =l= q(n,t);
xq1_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    xq1(n,t) =g= x(n,t) + q(n,t) - 1;
* -- product of x(n,t) and z(n,t) ---------------------------
zx_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     zx(n,t) =l= z(n,t); //t-1
zx_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     zx(n,t) =l= x(n,t)   ;
zx_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     zx(n,t) =g= x(n,t) + z(n,t) - 1;  // t-1
* -- product of q(n,t) and w1(n,t) --------------------------
w1q_1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    w1q(n,t) =l= M_ * (1-q(n,t));
w1q_2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    w1q(n,t) - w1(n,t) - M_ * (q(n,t)) =l= 0;
w1q_3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..    w1q(n,t) =g= wq(n,t);
* -- binary indicator x(n,t) when site n is infested first time in period t ----
x_u(n,t)$(t_ID(t)>1 and mask(n,t)>0)..      x(n,t) =g= (u(n,t)-u(n,t-1));
x_u1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     x(n,t) + u(n,t) + u(n,t-1) =l= 2;
x_u2(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     x(n,t) =l= u(n,t);
x_u3(n,t)$(t_ID(t)>1 and mask(n,t)>0)..     x(n,t) =l= 1-u(n,t-1);

* -- no treatments of sites n with the population densities < w_det in period t
q_v(n,t)$(t_ID(t)>1 and mask(n,t)>0)..      q(n,t) * w_det =l= w1(n,t);

* -- expectd pop.density in site n in period t vs t-1 without treatments -------
pop_nt0(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w1(n,t) =g= G(n) * w(n,t-1);

* -- pop.density in site n in period t vs. t-1 - old version, more non-zeroes
pop_nt(n,t)$(t_ID(t)>1 and mask(n,t)>0)..   w(n,t) =e=
 (G(n) * (wq(n,t) - wxq(n,t) + (1-e) * wq1(n,t)) + w_min * xq(n,t));
* -- pop.density in site n in period t vs. t-1 - new version - faster
pop_nt1(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =e=
  G(n) * (w(n,t-1) - wx(n,t) - wq1(n,t) * e + wxq1(n,t)) + w_min * (x(n,t)-xq1(n,t));
* -- pop.density in site n in period t vs. t-1 - new version, reduced
pop_nt1_(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =e=
  G(n) * (w(n,t-1) - wq1(n,t) * e) + w_min * (x(n,t)-xq1(n,t));

* -- pop.density in site n in period t vs. t-1 - spread-only scenario
pop_nt00(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =e=  G(n) * (w(n,t-1) - wx(n,t)) + w_min * x(n,t);
* -- pop.density in site n in period t vs. t-1 - spread-only scenario, new version, reduced
pop_nt00_(n,t)$(t_ID(t)>1 and mask(n,t)>0)..  w(n,t) =e=  G(n) * w(n,t-1) + w_min * x(n,t);

* -- spread from site n to other sites in period t, assign the pop.density= w_min to newly infested sites
buff_nodes(n,m,t)$(spread(n,m,t)>0 and t_ID(t)>1 and mask(n,t)>0).. u(m,t)$(spread(n,m,t)>0) =g= z(n,t);
buff_wmin(n,m,t)$(spread(n,m,t)>0 and t_ID(t)>1 and mask(n,t)>0)..  w(m,t)$(spread(n,m,t)>0) =g= zx(n,t) * w_min;
* -- x(n_t) = 0 if there is no new spread to n in period t from other sites
sum_buff_n(m,t)$(t_ID(t) > 1 and mask(m,t)>0).. x(m,t) =l= sum(n$(spread(n,m,t) > 0 and mask(n,t)>0), z(n,t));

* -- site n once infested remains infested
inf_u_t(n,t)$(t_ID(t)>1 and mask(n,t)>0).. u(n,t) =g= u(n,t-1);

* -- initialization to find feasible solution: a sequence of one-period solutions
* -- for periods 2,...,T:
* -- set the infested, treated status and pop.densities for periods t > t_ to 0
fix_q_t_0(n,t)$(t_ID(t) > t_step)..  q(n,t) =e= 0;

fix_u_t_(n,t)$((t_ID(t) < t_step  ) and t_ID(t)>1 ).. u(n,t) =e= u0(n,t);
fix_q_t_(n,t)$((t_ID(t) < t_step  ) and t_ID(t)>1 ).. q(n,t) =e= q0(n,t);
fix_w_t_(n,t)$((t_ID(t) < t_step  ) and t_ID(t)>1 ).. w(n,t) =e= w0(n,t);

* -- fix the infested, treated status and the pop.densities for periods 2,...,t-1
* -- so only one period t=card(t) need to be solved at a time
fix_u_t(n,t)$((ord(t) le (card(t)-1)) and t_ID(t)>1 ).. u(n,t) =e= u0(n,t);
fix_q_t(n,t)$((ord(t) le (card(t)-1)) and t_ID(t)>1 and q0(n,t)>0 ).. q(n,t) =e= q0(n,t);
fix_w_t(n,t)$((ord(t) le (card(t)-1)) and t_ID(t)>1 ).. w(n,t) =e= w0(n,t);
* -- enforcing the no-treatment spread only scenario
no_treat(n,t)$(mask(n,t)>0).. q(n,t) =e= 0;
* -- fix the treatment variable to q0(n,t)
fix_q(n,t)$(t_ID(t)>0 and mask(n,t)>0).. q(n,t) =e= q0(n,t);

* -- m_no = 2: full-size model - solved after solving the sequence of p_inf1 models
* -- min (detected inf. area * priority). s.t. budget constraints  over T periods
model p_inf   /
              obj_N_v
*              obj_N  fix_q
* test: calculating the budget based on pop.density just before the treatments
* ----- comment budget, wq2_1, wq2_2, wq2_3,   uncomment budgetG
              budgetG
              w_det1   w_det2   w_wspr   w_wspr1   w_wspr2   w_wmax   w_wmax1
              z_a_d3    pop_nt1     pop_nt0    q_v
              u_inf_nt  u_uninf_nt  w_inf_u    u_t       inf_u_t
              x_u       x_u1        x_u2       x_u3
              zx_1      zx_2        zx_3       xq1_1     xq1_2     xq1_3
              wx_1      wx_2        wx_3       wxq1_1    wxq1_2    wxq1_3
              wq1_1     wq1_2       wq1_3      wq_1      wq_2        wq_3
              buff_nodes       buff_wmin
             fix_q_t_0
              /
;
* -- m_no = 6: full-size model - solved after solving the sequence of p_inf1_cost models
* -- min (cost) + detected spread area penalty s.t. StS target
model p_inf_cost  /
              obj_cost  max_n_inf   V_budgetG
              w_det1    w_det2      w_wspr     w_wspr1   w_wspr2   w_wmax  w_wmax1
              z_a_d3    pop_nt1     pop_nt0    q_v
              u_inf_nt  u_uninf_nt  w_inf_u    u_t
              x_u       x_u1        x_u2       x_u3
              zx_1      zx_2        zx_3       xq1_1     xq1_2     xq1_3
              wx_1      wx_2        wx_3       wxq1_1    wxq1_2    wxq1_3
              wq1_1     wq1_2       wq1_3      wq_1      wq_2        wq_3
              buff_nodes       buff_wmin
              fix_q_t_0
              /
* One-period model to build a feasible solution - solve T one-period problems in sequence.
* Start from t_set=2 until t_set=t_ . After each solve starting from t=2,
* the model saves the decision variables, which are loaded on the next solving step, t=3,...,T
* After solving the initialization sequence for all t_ = T periods, the solutuion
* is used as a warm start the full problem p_inf_cost (m_no=6).
* -- m_no = 5: initialization - solve treatments for one period and spread for T periods
* -- min (cost) + detectd spread area penalty s.t. StS target
model p_inf1_cost  /
              obj_cost  max_n_inf   V_budgetG
              w_det1    w_det2    w_wspr   w_wspr1   w_wspr2  w_wmax    w_wmax1
              z_a_d3    pop_nt1     pop_nt0    q_v
              u_inf_nt  u_uninf_nt  w_inf_u    u_t
              x_u       x_u1        x_u2       x_u3
              zx_1      zx_2        zx_3       xq1_1     xq1_2     xq1_3
              wx_1      wx_2        wx_3       wxq1_1    wxq1_2    wxq1_3
              wq1_1     wq1_2       wq1_3      wq_1      wq_2        wq_3
              buff_nodes       buff_wmin
* -- t_set is set to t_ spread is solved for all peirods= t_set; the treatments are solved only for period tstep= card(t)
              fix_u_t_  fix_q_t_ fix_w_t_   fix_q_t_0
              /

* Obne-period  model to build a feasible solution - solves T one-period problems in sequence.
* Start from t_set=2 until t_set=t_ . After each solve starting from t=2,
* the model saves the decision variables, which are loaded on the next solving step, t=3,...,T
* After solving the initialization sequence for all t_ = T periods, the solutuion
* is used as a warm start the full problem p_inf (m_no=2).
* -- min (detected inf. area * priority) for one period t, tracking the spread for T periods
* -- m_no = 1: initialization - solve treatments for one period and spread for T periods
model p_inf1  /
              obj_N_v
              budgetG
              w_det1    w_det2      w_wspr     w_wspr1   w_wspr2   w_wmax   w_wmax1
              z_a_d3    pop_nt1     pop_nt0    q_v
              u_inf_nt  u_uninf_nt  w_inf_u    u_t       inf_u_t
              x_u       x_u1        x_u2       x_u3
              wq_1      wq_2        wq_3       wq1_1     wq1_2     wq1_3
              wxq1_1    wxq1_2      wxq1_3     wx_1      wx_2      wx_3
              xq1_1     xq1_2       xq1_3      zx_1      zx_2      zx_3
              buff_nodes            buff_wmin
* -- t_set is set to t_ spread is solved for all peirods= t_set; the treatments are solved for period tstep= card(t)
              fix_u_t_  fix_q_t_    fix_w_t_
              fix_q_t_0
              /

* Short-sighted problem when both treatments and spread are optimizaed for one period t only
* The one-period problem can be solved T periods in sequence starting from t_set=2 to t_set=t_ .
* After each solve starting from t=2, the model saves the decision variables and loads them
* on the next solving step, t=3,...,T, etc. After solving the initialization sequence
* for all t_ = T periods, the solutuion can be used as a warm start the full problem p_inf.
* -- min detected inf. area * priority for one period, spread and treatments for one period
* -- m_no = 4: initialization - solve treatments and spread for 1 period
model p_inf2  /
              obj_N_v1
              budgetG
              w_det1    w_det2      w_wspr     w_wspr1   w_wspr2   w_wmax   w_wmax1
              z_a_d3    pop_nt1     pop_nt0    q_v
              u_inf_nt  u_uninf_nt  w_inf_u    u_t       inf_u_t
              x_u       x_u1        x_u2       x_u3
              wq_1      wq_2        wq_3       wq1_1     wq1_2     wq1_3
              wxq1_1    wxq1_2      wxq1_3     wx_1      wx_2      wx_3
              xq1_1     xq1_2       xq1_3      zx_1      zx_2      zx_3
              buff_nodes            buff_wmin
* -- t_set is set to t_ spread is solved for all peirods= t_set; the treatments are solved for period tstep= card(t)
              fix_u_t_  fix_q_t_    fix_w_t_
              fix_q_t_0

              /

* -- m_no = 3: spread model only, no treatments - solve first before the treatment models
model p_inf_s /
              obj_N_v
              w_det1    w_det2
              w_wspr    w_wspr1     w_wspr2    w_wmax    w_wmax1
              z_a_d1    z_a_d2      z_a_d3     pop_nt00
              u_inf_nt  u_uninf_nt  w_inf_u    u_t0
              wx_1      wx_2        wx_3
              x_u       x_u1        x_u2       x_u3
              zx_1      zx_2        zx_3
              wq_1      wq_2        wq_3
              buff_nodes       buff_wmin       no_treat
              sum_buff_n
              /
;

Option MIP = GUROBI;
p_inf.optfile = 1;     p_inf1.optfile = 1;      p_inf2.optfile = 1;
p_inf_s.optfile = 1;   p_inf_cost.optfile = 1;  p_inf1_cost.optfile = 1;

* -- spread-only model - solve for t_set=tstep=card(t_) prior to solving the treatment model
if (m_no eq 3, solve p_inf_s minimize obj_val using mip; NZ=p_inf_s.numNZ/10e6;) // m_no = 3;

* -- one-period modle - solve p_inf2 in sequence for t_set=tstep=2 ... card(t_)
if (m_no eq 4, solve p_inf2 minimize obj_val using mip; NZ=p_inf2.numNZ/10e6; P_names=1; if (t_step>2, P_names=0;));

* -- initialization: solve p_inf1 for a sequence of tstep=2 ... tstep=card(t_)
if (m_no eq 1, solve p_inf1 minimize obj_val using mip; NZ=p_inf1.numNZ/10e6; P_names=1; if (t_step>2, P_names=0;));

* -- full model - solve p_inf
if (m_no eq 2, solve p_inf minimize obj_val using mip;  NZ=p_inf.numNZ/10e6;  );

* -- initialization: solve p_inf_cost to minimize the treatment cost for the StS target
if (m_no eq 5, solve p_inf1_cost minimize obj_val using mip; NZ=p_inf1_cost.numNZ/10e6; P_names=1; if (t_step>2, P_names=0;));// solving the sequential problem (first)

* -- full model with the StS target - solve after solving the sequence of p_inf1_cost
if (m_no eq 6, solve p_inf_cost minimize obj_val using mip; NZ=p_inf_cost.numNZ/10e6; P_names=1;);

Scalar telapsed;
Parameter gap_v, obj_e; // Retrieve the optimality gap

if (m_no eq 1, obj_e = abs(p_inf1.objest); gap_v = (abs((obj_e - p_inf1.objval) / obj_e)););
if (m_no eq 4, obj_e = abs(p_inf2.objest); gap_v = (abs((obj_e - p_inf2.objval) / obj_e)););
if (m_no eq 2, obj_e = abs(p_inf.objest);  gap_v = (abs((obj_e - p_inf.objval)  / obj_e)););
if (m_no eq 3, obj_e = abs(p_inf_s.objest);  gap_v = (abs((obj_e - p_inf_s.objval)  / obj_e)););
if (m_no eq 5, obj_e = abs(p_inf1_cost.objest); gap_v = (abs((obj_e - p_inf1_cost.objval) / obj_e)););
if (m_no eq 6, obj_e = abs(p_inf_cost.objest); gap_v = (abs((obj_e - p_inf_cost.objval) / obj_e)););

telapsed = TimeElapsed;

u0(n,t) = u.l(n,t);       w0(n,t) = w.l(n,t);       v0(n,t) = v.l(n,t);
q0(n,t) = q.l(n,t);       z0(n,t) = z.l(n,t);       x0(n,t) = x.l(n,t);
d0(n,t) = d.l(n,t);       a0(n,t) = a.l(n,t);
xq10(n,t) = xq1.l(n,t);   w10(n,t) = w1.l(n,t);
wx0(n,t) = wx.l(n,t);     zx0(n,t) = zx.l(n,t);     wxq10(n,t) = wxq1.l(n,t);
wq10(n,t) = wq1.l(n,t);   wq0(n,t) = wq.l(n,t);     wq20(n,t) = wq2.l(n,t);

* -- calculating the mask of eligible sites for spread-only scenario
* -- the mask is then used in initialization and treatment scenarios
file f2c / '%outpath%/mask_spread.txt' /;   f2c.pc = 6; f2c.nd = 0; f2c.pw=32767;
file fw0 / '%outpath%/_msk.txt' /;  fw0.pc = 6; fw0.nd = 0; fw0.pw=32767;

Parameter mask_tmp(n), mask_(n,t), spread_nm(n,m);

* -- writing the reduced mask based on spread-only scenario - sites u(n,t) with 20000-m buffer around
if (m_No eq 3,
    loop(n, loop(t, mask_(n,t)$(u0(n,t)>0) = u0(n,t);
    if(u0(n,t)>0, loop(m$(dist(n,m) < buffer_mask), mask_(m,t) = 1;)); ))
    loop(n, loop(t$(t_ID(t)>1), if(mask_(n,t-1) > 0, mask_(n,t) = 1;) ); );
    loop(n, loop(t, if(mask_(n,t)>0, put f2c n.tl t.tl mask_(n,t) /;)));
    loop(n,put fw0 n.tl; loop(t,put fw0 mask_(n,t)); put fw0 /;);
   );
putclose f2c;
putclose fw0;

scalar tset; tset=%t_set%;
* -- calculating the summaries -------------------------------------------------
Parameter n_treat(t), t_cost(t), n_inf(t), n_spr(t), n_det(t), w_tr_av(t), n_1st(t);
w_tr_av(t)$(sum(n, q0(n,t))>0) = sum(n, q0(n,t) * w10(n,t)) / sum(n, q0(n,t)); // mean pop.dens. of the treated sites in period t
n_treat(t) = sum(n, q0(n,t)); // number of treatead sites in period t
t_cost(t) = sum(n, q0(n,t) * (cost_fix + cost_pop * G(n) * wq10(n,t) )); // budget(n) calculations
n_inf(t) = sum(n, u0(n,t));
n_spr(t) = sum(n, z0(n,t));
n_det(t) = sum(n, v0(n,t));
n_1st(t) = sum(n, x0(n,t));

Parameter G_av, w_spr_av, x_av, x_, v_av;
G_av = sum(n, G(n)) / card(n);
w_spr_av = sum(n, w_spr(n)) / card(n);

if(card(t)>2, x_av = sum(t$(t_ID(t)>2), n_1st(t)) / (t_step - 2));
x_ = sum(t$(t_ID(t)>1 and t_ID(t) le t_step), n_1st(t)) / (t_step - 1);
if(card(t)>2, v_av = sum(t$(t_ID(t)>2 and t_ID(t) le t_step), sum(n, (v0(n,t) - v0(n,t-1)) )) / (t_step - 2););

Parameter t_inf(n), t_inf_det(n), t_delay(n), t_delay1(n), delay;
* --- period t when the infested site n becomes detectable
file f2a / '%outpath%/inf_t_det.txt' /;   f2a.pc = 6; f2a.nd = 0; f2a.pw=32767; // period t when infestation becomes detectable

Parameter det_year(n), tmp_sw, tmp_value;
loop(n, tmp_sw = 0; tmp_value = 0;
         if (w0(n,"1") ge w_det, tmp_value = t_ID("1"); tmp_sw=1;);
         loop(t$(t_ID(t)>1 and t_ID(t) le t_step),
         if(v0(n,t) eq 1 and v0(n,t-1) eq 0 and tmp_sw eq 0, tmp_sw = 1; tmp_value = t_ID(t); ); );
put f2a n.tl tmp_value /;  t_inf_det(n) = tmp_value;
                   );
putclose f2a;
* --- period t when the site n becomes is infested
file f2b / '%outpath%/inf_t.txt' /;   f2b.pc = 6; f2b.nd = 0; f2b.pw=32767; // period t when site n is infested
loop(n, tmp_sw = 0; tmp_value = 0;
         if (w0(n,"1") > 0, tmp_value = t_ID("1"); tmp_sw=1;);
         loop(t$(t_ID(t)>1 and t_ID(t) le t_step),
         if(x0(n,t) eq 1 and tmp_sw eq 0, tmp_sw = 1; tmp_value = t_ID(t); ); );

put f2b n.tl tmp_value /;  t_inf(n) = tmp_value;
                       );
putclose f2b;

t_delay1(n)$(t_inf_det(n) > t_inf(n) and t_inf_det(n) > 1) = 1;
t_delay(n)$(t_inf_det(n) > t_inf(n) and t_inf_det(n) > 1) = t_inf_det(n) - t_inf(n);
delay$(sum(n,t_delay1(n))>0) = sum(n,t_delay(n)) / sum(n,t_delay1(n));

Parameter retreat(n), treat_T(n), retreat_1, retreat_av, obj_v0, obj_noprior;
treat_T(n) = sum(t$(t_ID(t) le t_step), u0(n,t));
retreat(n)$(treat_T(n)>0) = sum(t$(t_ID(t) le t_step), q0(n,t)) / treat_T(n);
retreat_1 = sum(n$(retreat(n)>0), 1);
retreat_av$(retreat_1>0) = sum(n, retreat(n)) / retreat_1;

obj_v0= sum(n, prior(n) * v0(n,"2")) +
        sum(t$(t_ID(t)>2 and t_ID(t) le t_step ), sum(n, prior(n) * (v0(n,t) - v0(n,t-1)) ));
obj_noprior= sum(n, v0(n,"2")) +
        sum(t$(t_ID(t)>2 and t_ID(t) le t_step ), sum(n, (v0(n,t) - v0(n,t-1)) ));

Parameter n_inf_0, n_inf_, n_inf_s, StS_, B_av;
n_inf_ =  sum(n$(mask(n,"2")>0), prior(n) * v0(n,"2"))+
          sum(t$(t_ID(t) le t_step and t_ID(t)>2),
          sum(n$(mask(n,t)>0), prior(n) *  (v0(n,t) - v0(n,t-1))   ));

n_inf_s = sum(t$(t_ID(t) le t_step and t_ID(t)>1),max_obj_t(t) );
n_inf_0 = sum(t$(t_ID(t) le t_step and t_ID(t)>1),max_obj_t(t) ) * StS_target;
B_av = sum(t$(t_ID(t)>1 and t_ID(t) le t_step), b_.l(t)) / (t_step - 1);
StS_$(n_inf_s > 0) = n_inf_ / n_inf_s;

* -- wrting the objective values for periods t in the spread-only scenario
Parameter max_obj_t_s(t);
max_obj_t_s("2") = sum(n$(mask(n,"2")>0), prior(n) * v0(n,"2") ) ;
max_obj_t_s(t)$(t_ID(t)>2) = sum(n$(mask(n,t)>0), prior(n) * (v0(n,t) - v0(n,t-1)) );

file fwa / '%outpath%/max_obj_t.txt' /;  fwa.pc = 6; fwa.nd = 5; fwa.pw=32767;
if (m_no eq 3,    loop(t$(t_ID(t)>1), put fwa t.tl max_obj_t_s(t) /;   );  );
putclose fwa;
file fwb / '%outpath%/obj_t.txt' /;  fwb.pc = 6; fwb.nd = 5; fwb.pw=32767;
loop(t$(t_ID(t)>1), put fwb t.tl max_obj_t_s(t) /;   ); putclose fwb;

Scalar t0step;  t0step = %t_out%; // Number ot time steps to output in results1.txt
* -- writing the summary results -----------------------------------------------
file f_ / '%outpath%/results1.txt' /; f_.pc = 6; f_.nd = 3; f_.ap = 1; f_.pw=32767;

if(P_names>0,
put f_ "Stime" "gap" "s_time   " "obj_f" "objv0" "n_det" "M_NZ" "m_No" "N" "T" "Tstep";
put f_ "MaxD" "B" "StS" "n_inf" "n0inf" "e" "G" "x_" "x_t" "v_av" "costF" "costP";
put f_ "delay" "trFrq" "KernC" "w_min" "w_det" "w_spr" "w_max";
put f_ "w_tr:"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ "n_tr:"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ "cost:"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ "n_1st"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ "n_inf"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ "n_spr"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ "n_det"; loop(t$(t_ID(t) le t0step), put f_ t.tl);
put f_ /;
  );
f_.nd = 1; put f_ telapsed;      f_.nd = 3; put f_ gap_v ;
f_.nd = 0; put f_ system.time;
f_.nd = 2; put f_ obj_val.l obj_v0 obj_noprior NZ;
f_.nd = 0; put f_ m_no card(n) card(t) t_step Max_dist;   f_.nd = 2;
if (m_No eq 5 or m_No eq 6, put f_ B_av; else put f_ B;);
f_.nd = 3; put f_ StS_;
f_.nd = 2; put f_ n_inf_ n_inf_0 e G_av x_ x_av v_av cost_fix cost_pop;
f_.nd = 4; put f_ delay retreat_av;
if(gen_spread eq 1, f_.nd = 5; put f_ Kern_coef;);
if(gen_spread eq 0, put f_ "file";);
f_.nd = 4; put f_ w_min w_det w_spr_av w_max;
f_.nd = 4; put f_ "w_tr:"; loop(t$(t_ID(t) le t0step), put f_ w_tr_av(t));
f_.nd = 0; put f_ "n_tr:"; loop(t$(t_ID(t) le t0step), put f_ n_treat(t));
f_.nd = 3; put f_ "cost:"; loop(t$(t_ID(t) le t0step), put f_ t_cost(t));
f_.nd = 0; put f_ "n_1st"; loop(t$(t_ID(t) le t0step), put f_ n_1st(t));
f_.nd = 0; put f_ "n_inf"; loop(t$(t_ID(t) le t0step), put f_ n_inf(t));
f_.nd = 0; put f_ "n_spr"; loop(t$(t_ID(t) le t0step), put f_ n_spr(t));
f_.nd = 0; put f_ "n_det"; loop(t$(t_ID(t) le t0step), put f_ n_det(t));


put f_ /;  putclose f_;


* -- output file for pasting to Excel table - test_spread_10.xlsx and test_spread_15.xlsx
file f2 / '%outpath%/popul_nt.txt' /;   f2.pc = 6; f2.nd = 4; f2.pw=32767;
put f2 "n_ID" "w_min:" w_min "w_det" w_det "w_spr" w_spr_av "w_max"  w_max  ;
put f2 "Gmax" G0 "G_" G_av "x_t2" x_av "x" x_ "v" v_av "e" e "B" B "Max_d" Max_dist "M_" M_ /;
f2.nd = 0;

loop(n, put f2 n.tl;   loop(t$(t_ID(t) le t0step), put f2 u0(n,t));
put f2 "w:"; loop(t$(t_ID(t) le t0step), if(w0(n,t) eq 0, f2.nd = 0;); if(w0(n,t) > 0, f2.nd = 4;);
if (q0(n,t) eq 1, put f2 (w10(n,t)*(-1)););
if (q0(n,t) eq 0,put f2 w0(n,t)));         f2.nd = 0;
put f2 "w1(q):"; loop(t$(t_ID(t) le t0step), if(w10(n,t) eq 0, f2.nd = 0;); if(w10(n,t) > 0, f2.nd = 4;);
if (q0(n,t)>0, put f2 w10(n,t)); if (q0(n,t) eq 0, f2.nd = 0; put f2 q0(n,t));
); f2.nd = 0;

put f2 "v:"; loop(t$(t_ID(t) le t0step), put f2 v0(n,t));
put f2 "z:"; loop(t$(t_ID(t) le t0step), put f2 z0(n,t));
put f2 "x:"; loop(t$(t_ID(t) le t0step), put f2 x0(n,t));
put f2 "q:"; loop(t$(t_ID(t) le t0step), put f2 q0(n,t));

put f2 /; );   putclose f2;

$ontext
file f0 / '%outpath%/indicators_nt.txt' /;   f0.pc = 6; f0.nd = 0; f0.pw=32767;

loop(t,
put f0 "t" "n" " " "w_min" "w_det" "w_spr" "w_max" "W1nt" "Wnt" ;
if(t_ID(t)>1, put f0 "Wnt-1" "u_t-1"); if(t_ID(t) eq 1, put f0 "  " "   ");
put f0 "u" "v" "x" " " "a" "d" "z";
if(t_ID(t)>1, put f0 "__" "q_t-1" "Wnt-1" "WQnt" "WQ1nt" );
put f0 "  " "q" "x" "xq";
if(t_ID(t)>1, put f0 "__" "q" "x" "Wnt-1" "WXQnt" "WXnt" );
put f0 "x" "z" "zx" /;

loop(n,
put f0 t.tl n.tl; f0.nd=4; put f0 "w" w_min w_det w_spr w_max w1.l(n,t) w0(n,t);
if(t_ID(t)>1, put f0 w0(n,t-1) u.l(n,t-1)); if(t_ID(t) eq 1, put f0 "__" "   "); f0.nd=0;
put f0 u.l(n,t) v.l(n,t) x.l(n,t)    "z=ad"  a.l(n,t) d.l(n,t) z.l(n,t) ;
if(t_ID(t)>1, put f0 "wqwq1" q.l(n,t-1); f0.nd=4; put f0  w.l(n,t-1) wq.l(n,t) wq1.l(n,t)); f0.nd=0;
put f0 "xq wx" q.l(n,t) x.l(n,t) xq1.l(n,t);
if(t_ID(t)>1,
put f0 "wxq1" q.l(n,t) x.l(n,t); f0.nd=4; put f0 w.l(n,t-1) wxq1.l(n,t) wx.l(n,t) ;
 );
f0.nd=0; put f0  x.l(n,t) z.l(n,t) zx.l(n,t);

put f0 /;) put f0 /;);   putclose f0;
$offtext


scalar incr; incr = 0;
* -- pop.densities just before treatments ---------------
file f13 / '%outpath%/treat_dens.txt' /;   f13.pc = 6; f13.nd = 5; f13.pw=32767;
loop(t, loop(n, if(q0(n,t)>0, put f13 n.tl t.tl w10(n,t) /; incr=incr+1; )));
if(incr eq 0, loop(n$(ord(n) eq 1), put f13 n.tl "1" q0(n,"1") /;) );
incr = 0;
putclose f13;

* -- writing decision variables ----------
file f3 / '%outpath%/q.txt' /;   f3.pc = 6; f3.nd = 0; f3.pw=32767;
loop(n, loop(t, if(q0(n,t) > 0, put f3 n.tl t.tl q0(n,t) /; incr=incr+1; )));
if(incr eq 0, loop(n$(ord(n) eq 1), put f3 n.tl "1" q0(n,"1") /;) );
incr = 0;    putclose f3;

file f4 / '%outpath%/w.txt' /;   f4.pc = 6; f4.nd = 6; f4.pw=32767;
loop(n, loop(t, if(w0(n,t) > 0, put f4 n.tl t.tl w0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4 n.tl "1" w0(n,"1") /;) );
incr = 0; putclose f4;

file f4a / '%outpath%/w1.txt' /;   f4a.pc = 6; f4a.nd = 6; f4a.pw=32767;
loop(n, loop(t, if(w10(n,t) > 0, put f4a n.tl t.tl w10(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4a n.tl "1" w10(n,"1") /;) );
incr = 0; putclose f4a;

file f4b / '%outpath%/wx.txt' /;   f4b.pc = 6; f4b.nd = 6; f4b.pw=32767;
loop(n, loop(t, if(wx0(n,t) > 0, put f4b n.tl t.tl wx0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4b n.tl "1" wx0(n,"1") /;) );
incr = 0; putclose f4b;

file f4c / '%outpath%/wq.txt' /;   f4c.pc = 6; f4c.nd = 6; f4c.pw=32767;
loop(n, loop(t, if(wq0(n,t) > 0, put f4c n.tl t.tl wq0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4c n.tl "1" wq0(n,"1") /;) );
incr = 0; putclose f4c;

file f4d / '%outpath%/wq1.txt' /;   f4d.pc = 6; f4d.nd = 6; f4d.pw=32767;
loop(n, loop(t, if(wq10(n,t) > 0, put f4d n.tl t.tl wq10(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4d n.tl "1" wq10(n,"1") /;) );
incr = 0;  putclose f4d;

file f4e / '%outpath%/wq2.txt' /;   f4e.pc = 6; f4e.nd = 6; f4e.pw=32767;
loop(n, loop(t, if(wq20(n,t) > 0, put f4e n.tl t.tl wq20(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4e n.tl "1" wq20(n,"1") /;) );
incr = 0; putclose f4e;

file f4f / '%outpath%/wxq1.txt' /;   f4f.pc = 6; f4f.nd = 6; f4f.pw=32767;
loop(n, loop(t, if(wxq10(n,t) > 0, put f4f n.tl t.tl wxq10(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f4f n.tl "1" wxq10(n,"1") /;) );
incr = 0; putclose f4f;

file f5 / '%outpath%/z.txt' /;   f5.pc = 6; f5.nd = 0; f5.pw=32767;
loop(n, loop(t, if(z0(n,t) > 0, put f5 n.tl t.tl z0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f5 n.tl "1" z0(n,"1") /;) );
incr = 0; putclose f5;

file f6 / '%outpath%/x.txt' /;   f6.pc = 6; f6.nd = 0; f6.pw=32767;
loop(n, loop(t, if(x0(n,t) > 0, put f6 n.tl t.tl x0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f6 n.tl "1" x0(n,"1") /;) );
incr = 0; putclose f6;

file f7 / '%outpath%/u.txt' /;   f7.pc = 6; f7.nd = 0; f7.pw=32767;
loop(n, loop(t, if(u0(n,t) > 0, put f7 n.tl t.tl u0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f7 n.tl "1" u0(n,"1") /;) );
incr = 0; putclose f7;

file f8 / '%outpath%/v.txt' /;   f8.pc = 6; f8.nd = 0; f8.pw=32767;
loop(n, loop(t, if(v0(n,t) > 0, put f8 n.tl t.tl v0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f8 n.tl "1" v0(n,"1") /;) );
incr = 0; putclose f8;

file f9 / '%outpath%/d.txt' /;   f9.pc = 6; f9.nd = 0; f9.pw=32767;
loop(n, loop(t, if(d0(n,t) > 0, put f9 n.tl t.tl d0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f9 n.tl "1" d0(n,"1") /;) );
incr = 0; putclose f9;

file f10 / '%outpath%/a.txt' /;   f10.pc = 6; f10.nd = 0; f10.pw=32767;
loop(n, loop(t, if(a0(n,t) > 0, put f10 n.tl t.tl a0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f10 n.tl "1" a0(n,"1") /;) );
incr = 0; putclose f10;

file f11 / '%outpath%/xq1.txt' /;   f11.pc = 6; f11.nd = 0; f11.pw=32767;
loop(n, loop(t, if(xq10(n,t) > 0, put f11 n.tl t.tl xq10(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f11 n.tl "1" xq10(n,"1") /;) );
incr = 0; putclose f11;

file f12 / '%outpath%/zx.txt' /;   f12.pc = 6; f12.nd = 0; f12.pw=32767;
loop(n, loop(t, if(zx0(n,t) > 0, put f12 n.tl t.tl zx0(n,t) /; incr=incr+1;)));
if(incr eq 0, loop(n$(ord(n) eq 1), put f12 n.tl "1" zx0(n,"1") /;) );
incr = 0; putclose f12;

file f14 / '%outpath%/w_spr.txt' /;   f14.pc = 6; f14.nd = 5; f14.pw=32767;
loop(n, if(w_spr(n) > 0, put f14 n.tl w_spr(n) / ));
putclose f14;
file f15 / '%outpath%/G.txt' /;   f15.pc = 6; f15.nd = 5; f15.pw=32767;
loop(n, if(G(n) > 0, put f15 n.tl G(n) / ));
putclose f15;

file f16 / '%outpath%/treat_freq.txt' /;   f16.pc = 6; f16.nd = 4; f16.pw=32767;
loop(n, if(retreat(n) > 0, put f16 n.tl retreat(n) / ));
putclose f16;

