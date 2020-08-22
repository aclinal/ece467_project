$ontext
         Load Flow Analysis in Polar Coordinates

         This program optimizes a dummy variable and solves the nonlinear
         power flow equations in polar form.
$offtext

*Define buses
Set i / 1, 2, 3, 4, 5 /;
alias(i,j);
Set Slackbus(i) / 1 /;
Set PVbus(i) / 2, 4 /;
Set PQbus(i) / 3, 5 /;

*Define table headers
Set BusHead / Pg, Pd, Qg, Qd /;
Set LimitHeads / Pmin, Pmax, Qmin, Qmax, Vmin, Vmax, deltaSlack /;
Set CapsHead / C /;
Set loadCurtailment / LC /;

*Set baseMVA to be 100MVA and line thermal limits at 15pu (assume 5000MVA limit)
Scalar BaseMVA / 100 /;
Scalar ThermalLimit / 500 /;

Table Bus(i,BusHead)     'Bus data'
         Pg      Pd      Qg      Qd
*        MW      MW    MVAR    MVAR
1             162.5              45
2       250     195              50
3         0   227.5       0      65
4       300     234              60
5         0   188.5       0      55
;

Table r(i,j)             'Resistance data'
         1       2       3       4       5
1        0    0.02    0.08       0       0
2     0.02       0    0.06    0.06    0.04
3     0.08    0.06       0    0.01       0
4        0    0.06    0.01       0    0.08
5        0    0.04       0    0.08       0
;

Table x(i,j)             'Reactance data'
         1       2       3       4       5
1        0    0.06    0.24       0       0
2     0.06       0    0.18    0.18    0.12
3     0.24    0.18       0    0.03       0
4        0    0.18    0.03       0    0.24
5        0    0.12       0    0.24       0
;

Table yCharging(i,j)     'yCharging data'
         1       2       3       4       5
1        0    0.03   0.025       0       0
2     0.03       0    0.02    0.02   0.015
3    0.025    0.02       0    0.10       0
4        0    0.02    0.10       0   0.025
5        0   0.015       0   0.025       0
;

*Note:   For PV buses, voltage specified to float between 0.98pu and 1.02pu
*        However, for PQ buses, Vmin=0.8 and Vmax=1.2 specified to direct the
*        nonlinear solution towards a 'realistic' one
Table Limits(i,LimitHeads)
         Pmin    Pmax    Qmin    Qmax    Vmin    Vmax    deltaSlack
1          -9       9      -9       9       1       1             0
2        0.55    0.55    -2.5     1.5    0.98    1.02
3      -2.275  -2.275   -0.65   -0.65    0.98    1.02
4        0.66    0.66      -3     1.8    0.98    1.02
5      -1.885  -1.885   -0.55   -0.55    0.98    1.02
;


Table Caps(i,CapsHead)   'Capacitor values'
         C
1        0
2        0
3        1
4        0
5      0.8
;

Table loadCurt(i,loadCurtailment)   'Load curtailment amount'
        LC
1    0.075
2     0.09
3    0.105
4    0.108
5    0.087
;



parameters
         g0(i,j) 'conductances'
         b0(i,j) 'susceptances'
         G(i,j) 'G-bus/Conductance matrix'
         B(i,j) 'B-bus/Susceptance matrix'
         YBusMag(i,j) 'Y-Bus magnitude'
         YBusAngle(i,j) 'Y-Bus angle'
         BusPU(i,BusHead)  'Bus powers, in PU'
;

BusPU(i,BusHead)         = Bus(i,BusHead)/BaseMVA;

g0(i,j) $(r(i,j) and x(i,j) ne 0)   =r(i,j)/(r(i,j)*r(i,j)+x(i,j)*x(i,j));
g0(i,j) $(r(i,j) and x(i,j) eq 0)   =0;
b0(i,j) $(r(i,j) and x(i,j) ne 0)   =(-1)*x(i,j)/(r(i,j)*r(i,j)+x(i,j)*x(i,j));
b0(i,j) $(r(i,j) and x(i,j) eq 0)   =0;

G(i,j) $(not sameas(i,j))        =-g0(i,j);
G(i,i)                           =sum(j,g0(i,j));
B(i,j) $(not sameas(i,j))        =-b0(i,j);
B(i,i)                           =sum(j,b0(i,j))+0.5*sum(j,yCharging(i,j));

YBusMag(i,j)                     =sqrt(G(i,j)*G(i,j)+B(i,j)*B(i,j));
YBusAngle(i,j) $(G(i,j) ne 0)    =arctan(B(i,j)/G(i,j));
YBusAngle(i,j) $(G(i,j) eq 0)    =0;
YBusAngle(i,j) $(not sameas(i,j) and YBusAngle(i,j) ne 0)    =YBusAngle(i,j)+3.14159;


variables
         D               'Dummy cost variable'
         P(i)            'Bus power net injection'
         Q(i)            'Bus reactive power net injection'
         V(i)            'Bus voltages'
         delta(i)        'Bus angles'
         Pgen(i)         'Generator real power output'
         Qgen(i)         'Generator reactive power output'
         PLine(i,j)      'Real power flow from bus i to j'
         QLine(i,j)      'Reactive power flow from bus i to j'
         SLine(i,j)      'Apparent power flow from bus i to j'
         PLoss           'Real losses'
         QLoss           'Reactive losses'
;

P.up(i)=Limits(i,'Pmax')+loadCurt(i,'LC');
P.lo(i)=Limits(i,'Pmin')+loadCurt(i,'LC');
Q.up(i)=Limits(i,'Qmax')+Caps(i,'C');
Q.lo(i)=Limits(i,'Qmin')+Caps(i,'C');
V.up(i)=Limits(i,'Vmax');
V.lo(i)=Limits(i,'Vmin');
delta.up(Slackbus)=Limits(Slackbus,'deltaSlack');
delta.lo(Slackbus)=Limits(Slackbus,'deltaSlack');
SLine.up(i,j)=ThermalLimit;

equation
         DUMMY                           'dummy objective function'
         Pflow(i)                        'net bus real power injection'
         Qflow(i)                        'net bus reactive power injection'
         PLineFct(i,j)                   'Real line flow'
         QLineFct(i,j)                   'Reactive line flow'
         SLineFct(i,j)                   'Apparent power flow'
         PLossFct                        'Compute real power losses'
         QLossFct                        'Compute reactive power losses'
         PgenFct(i)                      'Determine slackbus real power output'
         QgenFct(i)                      'Determine slackbus reactive power output'
;

DUMMY ..                         D =e= 1;
Pflow(i) ..                      P(i) =e= SUM(j,V(i)*V(j)*YBusMag(i,j)*cos(YBusAngle(i,j)+delta(j)-delta(i)));
Qflow(i) ..                      Q(i) =e= SUM(j,-V(i)*V(j)*YBusMag(i,j)*sin(YBusAngle(i,j)+delta(j)-delta(i)));
PLineFct(i,j) ..                 PLine(i,j)=e=V(i)*V(i)*g0(i,j)-V(i)*V(j)*g0(i,j)*cos(delta(j)-delta(i));
QLineFct(i,j) ..                 QLine(i,j)=e=-V(i)*V(i)*(b0(i,j)+(1/2)*yCharging(i,j))+V(i)*V(j)*b0(i,j)*sin(delta(j)-delta(i));
SLineFct(i,j) ..                 SLine(i,j)=e=sqrt((PLine(i,j)*PLine(i,j))+(QLine(i,j)*QLine(i,j)));
PLossFct ..                      PLoss =e= sum((i,j),PLine(i,j));
QLossFct ..                      QLoss =e= sum((i,j),QLine(i,j));
PgenFct(i) ..                    Pgen(i) =e= BaseMVA*(P(i)+BusPU(i,'Pd'));
QgenFct(i) ..                    Qgen(i) =e= BaseMVA*(Q(i)+BusPU(i,'Qd'));

option nlp=minos;
MODEL LOADFLOW /ALL/;
SOLVE LOADFLOW USING NLP MAXIMIZING D;

display BusPU;
display YBusMag;
display YBusAngle;

display G;
display B;
