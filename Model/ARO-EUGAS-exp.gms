
$TITLE ARO_GasMod

$ontext
--------------------------------------------------------------------------
..........................................................................

Version: final (v6)
Date: 16.03.2021

Ready:
- Level 2&3 correctly merged via with DUAL for lvl3 problem
- Column-and-constraint algorithm for lvl 1 <-> 2&3
- Uncertainty set for demand & supply (integer steps)
- shedding (ensure that shed_cost << M)
- excel-based scalable data input
- parametrization with pan-European data
- scenario matrix
- monthly data input
- storages in primal and dual form
- allows for both PCI list investments and 20% expantion of each link.

..........................................................................
--------------------------------------------------------------------------
$offtext

$eolcom #

$set LOADXLS

$setglobal  version         v6
$set        datadir         data\
$Set        DataIn          %version%_data_exp
$set        resultdir       results\
$set        result          %version%_result

*-------------------------------------------------------------------------
*                          Uncertainty budget settings
*-------------------------------------------------------------------------

SCALARS
M           /1e4/
GAMMA_D     /0/       #from 0 to card(d)x5
GAMMA_S     /0/       #from 0 to card(g)
;

*-------------------------------------------------------------------------
*                         Sets, mappings and model data
*-------------------------------------------------------------------------

SETS
    n               nodes
    g               gas production facilities
    d               demands
    mapG(g,n)       mapping (gen to n)
    mapD(d,n)       mapping (dem to n)

    l               arcs
    ex(l)           existing arcs (subset)
    pros(l)         prospective arcs (subset)
    mapSL(n,l)
    mapRL(l,n)

    t               /t1*t12/
    v               /v1*v6/

*   GAMMA_S+0  | GAMMA_S+1
    i Dsens         /i1*i5/
*   IC-cost 5 steps: [0.2 .. 1]
    j Ssens         /j1*j1/
;

ALIAS (n,nn);
ALIAS (v,vv);

Parameter
        LDATA(l,*)
        DDATA(d,*)
        GDATA(g,*)
        D_M_LO(d,t)
        D_M_DELTA(d,t)

        StorageData(d,*)
        EminS(d,t)
;

*-------------------------------------------------------------------------
*                               Reading data
*-------------------------------------------------------------------------

$onecho > ImportData.txt
set=n               rng=sets!B2       rdim=1
set=g               rng=sets!D2       rdim=1
set=d               rng=sets!F2       rdim=1
set=mapG            rng=sets!H2       rdim=2
set=mapD            rng=sets!L2       rdim=2
set=l               rng=network!B2    rdim=1
set=ex              rng=network!D2    rdim=1
set=pros            rng=network!F2    rdim=1
set=mapSL           rng=network!H3    rdim=2
set=mapRL           rng=network!K3    rdim=2
par=LDATA           rng=network!N2    rdim=1 cdim=1
par=DDATA           rng=demand!A1     rdim=1 cdim=1
par=GDATA           rng=supply!A1     rdim=1 cdim=1
par=D_M_LO          rng=demand!F1     rdim=1 cdim=1
par=D_M_DELTA       rng=demand!F34    rdim=1 cdim=1
par=StorageData     rng=Storage!A2    rdim=1 cdim=1
par=EminS           rng=Storage!K2    rdim=1 cdim=1
$offecho

$onUNDF
$call GDXXRW I=%datadir%%DataIn%.xlsx O=%datadir%%DataIn%.gdx cmerge=1 @ImportData.txt
$gdxin %datadir%%DataIn%.gdx
$onUNDF

$LOAD n,g,d,mapG,mapD
$LOAD l,ex,pros,mapSL,mapRL
$LOAD LDATA,DDATA,GDATA,D_M_LO,D_M_DELTA
$LOAD StorageData, EminS
$gdxin

*execute_unload "%version%_check.gdx"
*$stop

*-------------------------------------------------------------------------
*                               Variables
*-------------------------------------------------------------------------
Variables
    #MASTER
    Z_M             system costs (master)
    #SUB
    Z_D             system costs (dual)
    phiS(d,t)
;

Positive Variables
    #MASTER
    ETA             aux var to reconstruct obj. function of the ARO problem
    PG_M(g,t,v)     gas production level
    PL_M(l,t,v)     gas flows
    PLS_M(d,t,v)    gas shedding
    ES(d,t,v)       storage level
    PSC(d,t,v)      storage charge
    PSD(d,t,v)      storage discarge

    #SUB
    Pdem(d,t)       [RO] realization of demand
    PE(g)           [RO] realization of supply
    lam(n,t)        dual var
    phiPG(g,t)      dual var
    phiPL(l,t)      dual var
    phiLS(d,t)      dual var
    muSC(d,t)       dual var storage
    muSD(d,t)       dual var storage
    muSL(d,t)       dual var storage
    muSU(d,t)       dual var storage
    z_lam(n,t)      aux continuous var to linearize z(d)*lam(n) term
    z_phiPG(g,t)    aux continuous var to linearize zg(g)*phiPG(g) term
;


Binary Variables
    #MASTER
    inv_M(l)
    #SUB
    z(d,t)
    zg(g)
;

parameters
    Tol         / 1 /
    LB          / -1e10 /
    UB          / 1e10 /
    itaux       / 0 /
    IB          / 1000 /
;

Parameters
    Pdem_fixed(d,t,v)   [1 <- 2&3] fixed realisation of demand
    PE_fixed(g,t,v)     [1 <- 2&3] fixed realisation of supply
    inv_iter_hist(l,v)
    inv_cost_master
    report_main(j,i,*,*)
    report_decomp(v,i,j,*,*)
    report_UBperN(j,i,n,t,*)
    report_UBperG(j,i,g,t)
    IC(l)
*   inv_track(l)
*   count;
    ;

IC(l) = LDATA(l,"exp_costs_1");

*-------------------------------------------------------------------------
*                        EQUATIONS
*-------------------------------------------------------------------------
Equations
#   Master
    MP_OBJ
    MP_marketclear(n,t,v)
    MP_PG(g,t,v)         
    MP_PLEX(l,t,v)       
    MP_PLPROS(l,t,v)     
    MP_LS(d,t,v)         
    MP_ETA(v)
    MP_SCCAP
    MP_SDCAP
    MP_SLEV
    MP_SLEV0
    MP_SminL
    MP_SmaxL
    MP_IBudget

#   Sub
    DUAL_OBJ
    DUAl_PG
    DUAL_PL
    DUAL_PLS
    EQ_DEM(d,t)
    EQ_UB_D
    EQ_SUP(g)
    EQ_UB_S
    DUAL_SD
    DUAL_SC
    DUAL_ES
    DUAL_ES12
    lin_1, lin_2, lin_3, lin_4,
    lin_5, lin_6, lin_7, lin_8
;

*-------------------------------------------------------------------------
*                               MASTER PROBLEM
*-------------------------------------------------------------------------

MP_OBJ.. Z_M  =E=    SUM(l$pros(l), IC(l)*inv_M(l)) + ETA;

MP_IBudget.. SUM(l$pros(l), IC(l)*inv_M(l)) =l= IB;

MP_marketclear(n,t,vv)$(ord(vv) lt (itaux+1))..
                            SUM(g$mapG(g,n),  PG_M(g,t,vv))
                        -   SUM(l$mapSL(n,l), PL_M(l,t,vv))
                        +   SUM(l$mapRL(l,n), PL_M(l,t,vv))
                        -   SUM(d$mapD(d,n),  Pdem_fixed(d,t,vv) - PLS_M(d,t,vv))
                        -   SUM(d$mapD(d,n),  PSC(d,t,vv) - PSD(d,t,vv))
                            =e= 0;

MP_PG(g,t,vv)$(ord(vv) lt (itaux+1))..      PE_fixed(g,t,vv) - PG_M(g,t,vv)   =g= 0;

MP_PLEX(ex,t,vv)$(ord(vv) lt (itaux+1))..       LDATA(ex,'tr_limit')/12               - PL_M(ex,t,vv)   =g= 0;
MP_PLPROS(pros,t,vv)$(ord(vv) lt (itaux+1))..   inv_M(pros)*LDATA(pros,'tr_limit')/12 - PL_M(pros,t,vv) =g= 0;

MP_LS(d,t,vv)$(ord(vv) lt (itaux+1))..        DDATA(d,"dem_y_lo")*D_M_LO(d,t) - PLS_M(d,t,vv)  =g= 0;

*** Storage
MP_SCCAP(d,t,vv)$(ord(vv) LT (ITAUX+1))..  StorageData(d,'PCmax') - PSC(d,t,vv)  =g= 0;
MP_SDCAP(d,t,vv)$(ord(vv) LT (ITAUX+1)).. StorageData(d,'PDmax') - PSD(d,t,vv)  =g= 0;

MP_SLEV(d,t,vv)$((ord(vv) LT (ITAUX+1)) AND (ord(t) GT 1) AND (StorageData(d, 'EFFD') GT 0))..
                                            ES(d,t,vv)
                                        -   ES(d,t-1,vv)
                                        -   PSC(d,t,vv)*StorageData(d,'EFFC')
                                        +   PSD(d,t,vv)/StorageData(d,'EFFD')
                                            =e= 0;

MP_SLEV0(d,vv)$((ord(vv) LT (ITAUX+1)))..      ES(d,'t1',vv)
                                        -   StorageData(d,'E0')
                                        -   PSC(d,'t1',vv)*StorageData(d,'EFFC')
                                        +   PSD(d,'t1',vv)/StorageData(d,'EFFD')
                                            =e= 0;

MP_SminL(d,t,vv)$(ord(vv) LT (ITAUX+1)).. ES(d,t,vv) - EminS(d,t) =g= 0;
MP_SmaxL(d,t,vv)$(ord(vv) LT (ITAUX+1)).. StorageData(d,'EMAX') - ES(d,t,vv) =g= 0;

***
MP_ETA(vv)$(ord(vv) lt (itaux+1))..  ETA =G= SUM((g,t), GDATA (g,"P_cost")*PG_M(g,t,vv))
                                            +   SUM((l,t), LDATA (L,"tr_cost")*PL_M(l,t,vv))
                                            +   SUM((d,t), DDATA (d,"LS_cost")*PLS_M(d,t,vv));

*-------------------------------------------------------------------------
*                               SUB PROBLEM
*-------------------------------------------------------------------------

DUAL_OBJ.. Z_D =E=  + SUM((d,n,t)$mapD(d,n), + DDATA(d,'dem_y_lo')*D_M_LO(d,t)*lam(n,t)
                                             + DDATA(d,'dem_y_lo')*D_M_LO(d,t)*D_M_DELTA(d,t)*z_lam(n,t))
                    + SUM((g,t),             - GDATA(g,'p_y_up')/12*phiPG(g,t) + GDATA(g,'p_delta')/12*z_phiPG(g,t))
                    + SUM((ex,t),            - LDATA(ex,"tr_limit")/12*phiPL(ex,t))
                    + SUM((d,t),             - DDATA(d,"dem_y_lo")*D_M_LO(d,t)*phiLS(d,t))
                    + SUM((d,t),  - muSC(d,t)*StorageData(d,'PCmax')
                                  - muSD(d,t)*StorageData(d,'PDmax')
                                  + muSL(d,t)*EminS(d,t)
                                  - muSU(d,t)*StorageData(d,'EMAX'))
                    - SUM(d,      + StorageData(d,'E0')*phiS(d,'t1'))
                    ;

DUAl_PG(g,t)..    + sum(n$mapG(g,n),lam(n,t)) - phiPG(g,t)                             =l= + GDATA(g,"P_cost");

DUAL_PL(ex,t)..   - sum(n$mapSL(n,ex),lam(n,t)) + sum(n$mapRL(ex,n),lam(n,t)) - phiPL(ex,t) =l= + LDATA(ex,"tr_cost");

DUAL_PLS(d,t)..   + sum(n$mapD(d,n),lam(n,t))  - phiLS(d,t)                            =l= + DDATA(d,"LS_cost");

# UB
EQ_DEM(d,t)..   Pdem(d,t) =E= DDATA(d,'dem_y_lo')*D_M_LO(d,t) + DDATA(d,'dem_y_lo')*D_M_LO(d,t)*D_M_DELTA(d,t)*z(d,t);
EQ_UB_D..       SUM((d,t), z(d,t)) =L= GAMMA_D;

EQ_SUP(g)..     PE(g) =E= GDATA(g,'p_y_up') - GDATA(g,'p_delta')*zg(g);
EQ_UB_S..       SUM((g), zg(g)) =L= GAMMA_S;

# Storage
DUAL_SD(d,t).. + sum(n$mapD(d,n),lam(n,t)) - muSD(d,t) =l= + 1/StorageData(d,'EFFD')*phiS(d,t);   #deriv PSD

DUAL_SC(d,t).. - sum(n$mapD(d,n),lam(n,t)) - muSC(d,t) =l= - StorageData(d,'EFFC')*phiS(d,t);     #deriv PSC

DUAL_ES(d,t)$(ord(t) LT 12).. - phiS(d,t) + phiS(d,t+1) + MUSL(d,t) - MUSU(d,t) =e= 0;            #deriv ES
DUAL_ES12(d,t)$(ord(t) EQ 12).. - phiS(d,t) + MUSL(d,t) - MUSU(d,t) =e= 0;                          #deriv ES

# Lin
lin_1(n,t).. -M*SUM(d$mapD(d,n), z(d,t))      =L=  z_lam(n,t)                           ;
lin_2(n,t)..  z_lam(n,t)                      =L=  M * SUM(d$mapD(d,n), z(d,t))         ;
lin_3(n,t).. -M*(1-SUM(d$mapD(d,n), z(d,t)))  =L=  lam(n,t) - z_lam(n,t)                ;
lin_4(n,t)..  lam(n,t) - z_lam(n,t)           =L=  M * (1 - SUM(d$mapD(d,n), z(d,t)))   ;

lin_5(g,t)..  -M*zg(g)                        =L=  z_phiPG(g,t)               ;
lin_6(g,t)..  z_phiPG(g,t)                    =L=  M * zg(g)                ;
lin_7(g,t)..  -M*(1-zg(g))                    =L=  phiPG(g,t) - z_phiPG(g,t)  ;
lin_8(g,t)..  phiPG(g,t) - z_phiPG(g,t)       =L=  M * (1 - zg(g))          ;

$include model_eq_stor.gms
*-------------------------------------------------------------------------
*                               Column-and-constraint decomposition
*-------------------------------------------------------------------------
* to allow CPLEX correctly detect rays in an infeasible problem
* optimality and feasibility tolerances are very small to avoid primal degeneration
*file COPT / cplex.opt /
*put COPT putclose 'ScaInd 0' / 'LPMethod 1' / 'PreInd 0' / 'EpOpt 1e-9' / 'EpRHS 1e-9' / ;
*ARO_MASTER.OptFile = 1;
*ARO_LVL23.OptFile = 1;

OPTION optcr=0;
OPTION MIP = CPLEX;
OPTION ResLim = 36000;

* option solveopt = clear;
ARO_MASTER.solveOpt = 2;
ARO_LVL23.solveOpt = 2;

*--------------------------------------------
#Settings of the test run
*itaux  = 1;
*Pdem_fixed(d,t,v) = DDATA(d,'dem_y_lo')*D_M_LO(d,t) ;
*PE_fixed(g,t,v)   = GDATA(g,'p_y_up')/12;
*inv_M.fx(l) = 0;

*SOLVE ARO_MASTER_V0 USING MIP MINIMIZING Z_M;
*execute_unload "%resultdir%%version%_masterV0.gdx";

*ARO_MASTER.savepoint=1;
*SOLVE ARO_MASTER USING MIP MINIMIZING Z_M;
*execute_unload "%resultdir%%version%_master.gdx";

*SOLVE ARO_LVL23 USING MIP MAXIMIZING Z_D;
*execute_unload "%resultdir%%version%_sub.gdx";

*$stop

*--------------------------------------------

loop(j,

    GAMMA_D = 0;
    GAMMA_S = 0;
    inv_iter_hist(l,v)  = 0;
    LB                  = -1e10;
    UB                  =  1e10;
    itaux               = 0;

#IC(L) = (2*ord(j)/10)*LDATA(l,"exp_costs");
         IC(L) $ (ord(j) = 1) = LDATA(l,"exp_costs_3");
         IC(L) $ (ord(j) = 2) = LDATA(l,"exp_costs_1");
*         IC(L) $ (ord(j) = 3) = LDATA(l,"exp_costs_1");

    loop(i,

    inv_iter_hist(l,v)  = 0;
    LB                  = -1e10;
    UB                  =  1e10;
    itaux               = 0;
    inv_cost_master     = 0;

* ############################################################
* D_EXP scenario [2 rows (cost factos) x 4 columns]
$ontext
        GAMMA_D $ (ord(i) = 1) = 0;
        GAMMA_S $ (ord(i) = 1) = 0;
        GAMMA_D $ (ord(i) = 2) = 10;
        GAMMA_S $ (ord(i) = 2) = 0;
        GAMMA_D $ (ord(i) = 3) = 30;
        GAMMA_S $ (ord(i) = 3) = 0;
        GAMMA_D $ (ord(i) = 4) = 60;
        GAMMA_S $ (ord(i) = 4) = 0;
$offtext
* ############################################################
* D_EXP_IB scenario [1 row (cost factor) x 5 columns (IB)]

        GAMMA_D $ (ord(i) = 1) = 60;
        GAMMA_S $ (ord(i) = 1) = 0;
        IB $ (ord(i) = 1) = 20;

        GAMMA_D $ (ord(i) = 2) = 60;
        GAMMA_S $ (ord(i) = 2) = 0;
        IB $ (ord(i) = 2) = 40;

        GAMMA_D $ (ord(i) = 3) = 60;
        GAMMA_S $ (ord(i) = 3) = 0;
        IB $ (ord(i) = 3) = 60;

        GAMMA_D $ (ord(i) = 4) = 60;
        GAMMA_S $ (ord(i) = 4) = 0;
        IB $ (ord(i) = 4) = 80;

        GAMMA_D $ (ord(i) = 5) = 60;
        GAMMA_S $ (ord(i) = 5) = 0;
        IB $ (ord(i) = 5) = 100;

* ############################################################         

*    inv_track(l)=1;
*    count=0;

       loop (v $((UB - LB) GT Tol),     # stop if UB-LB is lower than Tol value
*      loop (v,                         # stop if vector of investments repeats

            Pdem_fixed(d,t,v) = DDATA(d,'dem_y_lo')*D_M_LO(d,t);
            PE_fixed(g,t,v)   = GDATA(g,'p_y_up')/12;

            itaux = ord(v);

*            if (ord(v)=1,
*                SOLVE ARO_MASTER_V0 USING MIP MINIMIZING Z_M;
*            else
                SOLVE ARO_MASTER USING MIP MINIMIZING Z_M;
*            );

            LB = Z_M.l;

            inv_iter_hist(l,v) = inv_M.l(l);
            inv_cost_master = SUM(l, IC(l)*inv_M.l(l));

            report_decomp(v,i,j,'itaux','')     = itaux;
            report_decomp(v,i,j,'GAMMA_D','')   = GAMMA_D                                       + EPS;
            report_decomp(v,i,j,'GAMMA_S','')   = GAMMA_S                                       + EPS;
            report_decomp(v,i,j,'IC-mult','')   = IC('l108')/LDATA('l108',"exp_costs_1")          + EPS;
            report_decomp(v,i,j,'PCI built',l)   = inV_M.l(l)                                        ;
            report_decomp(v,i,j,'PCI MEUR/y','') = inv_cost_master                              + EPS;
            report_decomp(v,i,j,'M_obj','')     = Z_M.L                                         + EPS;
            report_decomp(v,i,j,'ETA','')       = ETA.l                                         + EPS;
            report_decomp(v,i,j,'LB','')        = LB;
            report_decomp(v,i,j,'M_Shed','')    = SUM((d,t), PLS_M.l(d,t,v))                        + EPS;
            report_decomp(v,i,j,'M_GC','')      = SUM((g,t), GDATA (g,"P_cost")*PG_M.l(g,t,v))      + EPS;
            report_decomp(v,i,j,'M_LC','')      = SUM((l,t), LDATA (L,"tr_cost")*PL_M.l(l,t,v))     + EPS;
            report_decomp(v,i,j,'M_SC','')      = SUM((d,t), DDATA (d,"LS_cost")*PLS_M.l(d,t,v))    + EPS;

*           inv_track(l) = inv_M.l(l);
*           loop(l, if (inv_track(l) ne inv_M.l(l), count = count + 1); );
*           break$(count=0);
*           count=0;

$include    network_upd_exp.gms

            SOLVE ARO_LVL23 USING MIP MAXIMIZING Z_D;

            Pdem_fixed(d,t,v)    = Pdem.l(d,t);
            PE_fixed(g,t,v)      = PE.l(g)/12;

            UB = min([inv_cost_master + Z_D.l], UB);

            report_decomp(v,i,j,'network','')   = card(ex);
            report_decomp(v,i,j,'Sub_obj','')   = Z_D.L                       + EPS;
            report_decomp(v,i,j,'UB','')        = UB;
            report_decomp(v,i,j,'Sub_Shed','')  = SUM((d,t), DUAL_PLS.m(d,t)) + EPS;
            report_decomp(v,i,j,'UB-LB','')     = UB - LB                     + EPS;

*            execute_unload "%resultdir%ARO_%version%_sub_it.gdx";

$include    network_clean_exp.gms
        );

    report_main(j,i,'itaux','')   = itaux;
    report_main(j,i,'GAMMA_D','') = GAMMA_D              + EPS;
    report_main(j,i,'GAMMA_S','') = GAMMA_S              + EPS;
    report_main(j,i,'IC-mult','') = IC('l108')/LDATA('l108',"exp_costs_1") + EPS;
    report_main(j,i,'PCI built',l) = inV_M.l(l);
    report_main(j,i,'PCI MEUR/y','') = inv_cost_master    + EPS;
    report_main(j,i,'M_obj','')    = Z_M.L                + EPS;
    report_main(j,i,'ETA','')      = ETA.l                + EPS;
    report_main(j,i,'LB','')       = LB                   + EPS;
    report_main(j,i,'S_obj','')    = Z_D.L                + EPS;
    report_main(j,i,'UB','')       = UB                   + EPS;
    report_main(j,i,'UB-LB','')    = UB - LB              + EPS;

    report_UBperN(j,i,n,t,'Demand-peaks [per node]')     = lin_3.m(n,t) ;
    report_UBperN(j,i,n,t,'Supply-drops [sum per node]') = sum(g$mapG(g,n), lin_7.m(g,t))  ;
    report_UBperG(j,i,g,t)                               = lin_7.m(g,t) ;

* Standard steps D=10, S=1:
*GAMMA_D=GAMMA_D+10;
*GAMMA_S=GAMMA_S+1;
    );

*put_utility 'gdxout' / 'J-steps' j.tl:0;
);

*-------------------------------------------------------------------------
*                               Model output
*-------------------------------------------------------------------------

execute_unload "%resultdir%%version%_report(D-EXP-IB-final).gdx" report_main, report_decomp, report_UBperN, report_UBperG;
*execute_unload "%resultdir%%version%_solution(DS-exp).gdx"
