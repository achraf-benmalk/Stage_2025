Effect of Punching on the Lifetime of Polymer Pipes for
the Conveyance of Drinking Water

Mechanics of mAterials for enGineering and Integrity of Structures - MAGIS Paris

2023-2024

Seyedeh Saeideh SAADAT

Seyedeh.saadat@ens-paris-saclay.fr

Reviewer:

Cristian OVALLE RODAS

Supervisors:

Xavier COLIN

Juan Pablo MARQUEZ COSTA

Laboratoire Procédés et Ingénierie en Mécanique et Matériaux (PIMM)

51 Boulevard de l'Hôpital – 75013, Paris, France

2

Table of Content

1.

2.

2.1.

2.2.

Motivation ............................................................................................................................................................. 4

Bibliography .......................................................................................................................................................... 4
Polyethylene Pipe Properties ................................................................................................................................. 4

Behavior of PE Pipes Under Different Mechanical Loads .................................................................................... 6
Hydrostatic Internal Pressure Load ....................................................................................................................... 6

Failure Mechanisms Under Hydrostatic Pressure.................................................................................................. 6
Point Load ............................................................................................................................................................. 7

Point Load Test ..................................................................................................................................................... 8
Failure Mechanisms Under Point Loads ............................................................................................................... 8

Effect of Various Parameters on The Point Load on The Lifetime of HDPE Pipes .............................................. 9
Chemical Degradation Assessment: Effects of Chlorine Disinfection .................................................................10

2.3.

Antioxidant Concentration Profiles ......................................................................................................................10
Changes in Molecular Weight ..............................................................................................................................11

Crystallinity Profiles ............................................................................................................................................11
Influence on Mechanical Behavior .......................................................................................................................12

Fractography.........................................................................................................................................................12
Lifetime Prediction: Combined Influence Of Punching And Chemical Degradation ..........................................12

2.4.

3.

Scientific Approach ..............................................................................................................................................13
3.1.  Main ideas ............................................................................................................................................................13

3.2.

the Tools and the techniques ................................................................................................................................14
Phase 1: Computational Modelling ......................................................................................................................14

Phase 2: Chemical Degradation Assessment: Effects of Chlorine Disinfection ...................................................16
Results and discussion ..........................................................................................................................................19

4.

4.1.
4.2.

Phase 1: Computational Modelling ......................................................................................................................19
Phase 2: Chemical Degradation Assessment: Effects of Chlorine Disinfection ...................................................27

5.

Concluding Remarks ............................................................................................................................................31

3

NOMENCLATURE

ASTM
TGA
CRB
DSC
FEM
FEA
FNCT
FTIR
GPC
HDPE
HOCl
LEFM
LDPE
MRS
MW
MWD
NHP
NHR
NPT
OD
OIT
PE
PE100-RC
PENT
PLs
PLT
PLT+
SCC
SDR
SCG

American Society for Testing of Materials
Thermogravimetric analysis
Cracked Round Bar Test
Differential Scanning Calorimetry
Finite Element Method
Finite Element Analysis
Full Notch Creep Test
Fourier Transform InfraRed spectroscopy
Gel Permeation Chromatography
High-Density Polyethylene
hypochlorous acid
Linear Elastic Fracture Mechanics
Low-Density Polyethylene
Minimum Required Strength
Molecular Weight
Molecular Weight Distribution
Notched HDPE Pipes
Notched HDPE Ring
Notched Pipe Test
Outer Diameter
Oxidation Induction Time
Polyethylene
HDPE Resistant to Crack
Pennsylvania Notch Test
Point Load
Point Load Test
Accelerated Point Load Test
Stress Corrosion Cracking
Standard Dimension Ratio
Slow Crack Growth

4

1. Motivation
Among  the polymer materials used for fluid transport, High-Density  Polyethylene  has constituted a significant  portion of
conduits in drinking water distribution networks since the 1990s. HDPE is the most common pipe material utilized across
the world, especially in the production of piping materials. For over half a century this material has been choose for use in
drinking water systems due not only to its unique properties, like corrosion resistance, flexibility, light weight,  toughness,
easy installation, water neutrality and leak tightness, but also to its sufficiently  long service lifetime. Since these materials
were initially introduced, there has been a substantial improvement in their performance. As far as the durability of polymer
pipe and resistance to chlorine for the newest type of HDPE is estimated to be up to 100 years.

Nevertheless,  this  expectation  of  the  lifetime  of  polymer  pipes,  especially  PE pipes, is  assured  only  under  conditions  of
internal  overpressure  with  pure  distilled  water  (without  chlorinated  disinfectants)  under  pressure.  However,  we  should
consider there are other factors, which can reduce lifetime of the pipes.

Thanks to the new generation of HDPE, the conduit made of the last generation of resin has been known for its properties of
resistance to PL, PE100-RC, can be installed with trenchless techniques and coarse-grained backfill materials. Consequently,
there is a potential risk of the existence of a local inhomogeneous stress field caused by a rock pressing onto the outside wall
of a polymer pipe. Therefore, the service life of buried pipes can be significantly  reduced due to this punctual load during
non-conventional pipe installations. Therefore, with this stress concentration, we regard a premature failure in polymer pipes.

Otherwise, microscopic analyses carried out in various studies and articles have revealed that cracks caused by the effect of
PL appear on the internal surface of polymer conduits. This observation raises the question of whether the punching effect
is the sole cause of this damage, or whether  there is a complex interaction between the mechanical effect of punching and
the chemical degradation induced by chlorinated disinfectants, which are commonly used in the treatment of drinking water.

Guaranteeing the reliability of drinking water distribution networks is a significant challenge for public health and the proper
functioning  of urban infrastructure. Otherwise, every year, large investments  are buried in the ground when repairing and
extending drinking water networks The aim of our study is to precisely investigate the specific effect of puncturing on the
structural properties of polymer pipes used in drinking water systems. This study will help determine to what extent these
internal deformations affect the lifetime of the pipes and performance, and assess their potential to lead to system failures.
Additionally, it will explore whether the observed damage results solely from puncturing or is influenced by other factors,
such as chemical  degradation caused by chlorinated disinfectants.  Understanding  the combined impact of these factors is
crucial to better comprehend the mechanisms of degradation and to develop strategies for prevention and improvement.

Although  experimental  tests accounting  for these  PL effects  are rare due to their complexity  and  duration, therefore this
study focuses on the numerical prediction of PL effects over pipe life and their coupling with oxidative degradation initiated
by attack from disinfectants.

2. Bibliography

Introduction

In this section, we cover four subsections. We begin by introducing an overview of HDPE pipe properties, especially with
the  introduction  of  the  latest  HDPE  generation,  PE100.  Next,  we  explore  the  impact  of  different  loads  on  these  pipes,
especially PLs caused by rack trenchless installation. Additionally, we investigate the chemical degradation of polymer pipes
exposed  to  chlorinated  disinfectants  commonly  used  in  drinking  water  treatment.  Finally,  the  section  concludes  with  an
examination of lifetime prediction, specifically focusing on the combined influence of punching and chemical degradation
in polymer pipes.

2.1.  POLYETHYLENE PIPE PROPERTIES

PE piping materials consist of a polyethylene polymer with additives like colorants, stabilizers, and antioxidants to enhance
properties and protect during manufacturing, storage, and service. Classified as a thermoplastic, PE pipes can be fabricated
using heat and pressure, and field joints can be made through thermal melting processes, creating a permanent bond.

PE is also a semi-crystalline polymer, exhibiting ordered molecular structures with substantial portions aligning closely in
crystalline regions. This semi-crystalline nature results in a low glass transition temperature (Tg) of approximately -183°K (-
90°C), contributing to greater toughness, impact resistance, and resistance to rapid crack propagation [Handbook of PE Pipe,
2008]. The higher degree of crystallinity of a PE the higher the density and the better the chemical resistance. [ H. Bromstrup,
2004]

5

Different Generations of PE Pipes

In the first generation of PE production, the  objective was to create long chains, making the material vulnerable to point
loading and crack expansion due to the parallel arrangement of the chains. Sensitivity to notches and point loads led to pipe
failures, resolved only with the addition of comonomers in second generation PE. These comonomers branched the molecular
chains, significantly improving the resistance to crack expansion but still requiring costly sand coating due to the effects of
point  loading.  In  the  third  generation  of  PE,  unlike  second-generation,  comonomers  are  directly  integrated  into  the  long
chains,  introducing  branching  in  a  more  controlled  manner.  This  branching  results  in  substantially  enhanced  protection
against punctual loads [J. Lenz, 2001].

HDPE  evolves  from  the  first  generation  PE32,  PE40,  through  the  second  generation  PE63,  PE80,  the  third  generation
bimodal  PE100,  to  the  latest state-of-the-art  PE100  RC.  The  pipe  in  PE100-RC  has  the  same  characteristics  in  terms  of
density and weldability as an ordinary pipe in PE100 resin. They have also the same standard dimensions. The difference is
in the material's ability to resist of PL and slow cracking, which is almost 10 times greater and excellent mechanical properties
and environmental stress crack resistance [Z. Xiong, R. Wei, 2012]. The resistance to crack initiation in PE is associated with
its amorphous structure, which includes additional ramifications as shown in Figure.1.[Dossier technique of SUEZ]

PE100 Properties

Figure.1. Amorphous structure in different resins of polyethylene

PE100,  an  advanced  version  compared  to  the  previous  generation  of  PE,  is  predominantly  manufactured  using  HDPE,
resulting in higher density due to increased crystallization. Crystalline  PE, with a density of 1.0 g/cm³, contributes to the
material's strength, while the amorphous part has a density of 0.85 g/cm³. PE100's enhanced strength is achieved through a
higher rate of crystallization, creating a more rigid material with increased density.

The PE100-RC grade of PE, which is a third-generation PE pipe, has been designed, since the end of the twentieth century,
to have a bimodal MWD to achieve an optimal balance between strength, ductility, and environmental stress crack resistance.
This approach involves a mixture of two molecular sizes, resulting in improved  PL resistance and enhanced flexibility  in
finished  pipes compared to traditional single  modular PE100 material.  The first  molecule,  similar  to monomolecular  PE,
provides high  density,  stiffness,  and a corresponding high  hoop stress value, while  the second, more branched molecule,
improves notch stress resistance and flexibility. [ H. Bromstrup, 2004]

A comparative investigation on two bimodal pipe-grade PE with different PL resistance by Y. Huang, et al., indicates that
PE100-RC pipe with better SCG resistance had a stronger deformation recovery capability during stretching [Y. Huang et
al., 2020]. Moreover, the Effect of short-chain branches distribution on the fracture behavior of PE pipe resins by [X. He, et
al,. 2018] suggest that the chain structure, including MW and its distribution, plays a crucial role in the fracture behavior of
PE100-RC [X. He et al. 2018]

Poisson ratio for PE has been found to vary somewhat depending on the ultimate strain that is achieved, on temperature and
the density of the base resin. however, for typical working  stress, strain and temperature an value  between   0.3 to 0.45 is
applicable to all PE pipe materials regardless of their densities and also for both short and long-term in-service [ Handbook
of PE Pipe, 2008].

Toughness, ductility, resilience, and resistance to damage from external loads, vibrations, and pressure surges like water
hammer  are the  inherent  properties of the PE pipe. Even in cold weather,  they  are flexible  and can be handled  and bent
without issue. Based on results from internal pressure tests, the standard extrapolation method outlined in  [EN ISO 9080,
2003] categorizes these pipes based on their MRS, ensuring a minimum service life of 50 years. PE 100 is the most recently
developed PE grade and PE pipes with higher strength and toughness than earlier generation materials. PE 100 has an MRS
of 10.0MPa. [A. Frank et al.2009] To gather information on the resistance to these failure mechanisms, modern methods of
LEFM are used. [A. Frank et al.2009]

6
Young's modulus  and yield stress of HDPE can vary based on the supplier and the crystallinity  rate of the material. An
increase in the crystallinity rate generally leads to an increase in both the elasticity modulus and the yield stress. For PE100,
the Young's modulus ranges between 850-1700 MPa, and the yield stress varies between 22-26 MPa according to data sheets.

linear  coefficient of  thermal  expansion  for  polyethylene  varies  widely  according  to  different  sources.  For  low-density
polyethylene, it is approximately 2 x 10⁻⁴ °C⁻¹ at 20°C, but this coefficient increases to about 3.5 x 10⁻⁴ °C⁻¹ at 80°C. This
property is related to the specific volume  expansion,  which  for an isotropic material  is three times  the linear coefficient.
Another reference states a value of 18 x 10⁻⁴ °C⁻¹. In various technical data sheets from different suppliers, lower values as
low as 1.3 x 10⁻⁴ K⁻¹ can also be found.

2.2.  BEHAVIOR OF PE PIPES UNDER DIFFERENT MECHANICAL LOADS

The evaluation of PE pipe performance under varying mechanical loads is essential to guarantee the strength and
dependability of these piping systems. This section is dedicated to examining how PE pipes respond to two main types of
mechanical loads: pressure load and PL.

Hydrostatic Internal Pressure Load

Under  static  internal  pressure,  the  pipe  wall  experiences  a  triaxial  stress  state,  where  the  most  significant  stress  is  the
circumferential stress. The axial stress on the inner surface is only half as much as the circumferential stress, while the radial
stress corresponds to the internal pressure in the inner surface and on the outer surface, the radial stress is zero.

According to a standard elastic mechanics model, the radial stress (𝜎𝑟), circumferential stress (𝜎𝑡) and axial stress (𝜎𝑎) of a
pipe with inner pressure can be calculated according to these formulas:

𝜎𝑟 =

𝜎𝑡 =

2
𝑃𝑖 𝑟𝑖
2 −   𝑟𝑖
𝑟𝑜
2
𝑃𝑖 𝑟𝑖
2 −   𝑟𝑖
𝑟𝑜

2  (1 −

2  (1 +

𝜎𝑎 =

2
𝑃𝑖 𝑟𝑖
2
2 −   𝑟𝑖
𝑟𝑜

2
𝑟𝑜
𝑟2)

2
𝑟𝑜
𝑟2)

Eq.1.

Eq.2.

Eq.3.

Where  𝑟 is  the  radial  coordinate(mm),  𝑟𝑜  is  the  outer  radius(mm),  𝑟𝑖  is  the  inner  radius(mm)  and  𝑃𝑖  is  the  inner
pressure (N).

In  designing  pressure  pipes,  we  follow  the  normal  stress  hypothesis,  focusing  on  the  highest  normal  stress,
which  leads  to  hoop  stress.  When  a  pipe  is  exposed  to  internal  pressure  and  the  ratio  of  inside  to  outside
diameter ≤1.3, we specifically consider circumferential stress as:

𝜎ℎ𝑜𝑜𝑝 =

𝑃𝑖 𝐷𝑚
2𝑒

Eq.4.

Where 𝐷𝑚 is the average pipe diameter and 𝑒 is wall thickness (mm). If 𝐷𝑚 is replaced by 𝐷𝑜 − 𝑒 dimension formula is:

𝜎ℎ𝑜𝑜𝑝 =

𝑃𝑖
10

 .

 𝐷𝑜 − 𝑒
2𝑒

Eq.5.

Where 𝜎ℎ𝑜𝑜𝑝 is hoop stress or circumferential  stress (in N/mm2 or MPa), P = 𝑃𝑖 is the inner pressure (bar) and 𝐷𝑜is outer
diameter of the pipe (mm).

Failure Mechanisms Under Hydrostatic Pressure

The crucial quality of a plastic pressure pipe is its ability to withstand hydrostatic pressure, determining its durability under
internal pressure. As mentioned in the previous section, the three-dimensional stress in the pipe wall due to internal pressure
is  significant.  The  fracture  stress  in  plastic  pipes  under  internal  pressure  is  influenced  by  test  duration  and  temperature,
extensively studied since 1956. The results are depicted in a log-log graph of test stress against failure time.

The failure mechanisms  observed in internally  pressurized pipes have been illustrated in the stress-life curve of Figure.2.
The relationship between applied load and structural lifetime, reveals three distinct kinetic regimes [A. Frank et al., 2007].
Each of these three parts of the curve corresponds to a different failure mechanism.

7
In the first regime (A), at high stress levels and at relatively short times 𝑡𝑓, the material shows a ductile failure with large-
scale plastic deformation  and a gradual slope, which  depends upon the material's  density  [H. Bromstrup, 2004]; Usually,
plastic pipe systems are designed to operate below this region.

In the transition  to the second regime  (B), the steep part, the intermediate  stress range, the  curve  steepens, shifting  from
ductile to quassi-brittle failure. In this regime, the failure is characterized by creep crack initiation, creep crack growth, and
only small-scale crack tip plasticity. There is a widely accepted understanding that this failure region significantly impacts
the lifetime of applications intended for long-term  use [N. Brown and X. Lu, 1991,1993]. This knee-point - the transition
from the flat to the steep part of the curve - can only be seen for PE 100 materials, if at all, at elevated temperatures and very
long testing times (at 80 °C not before 10,000 hours). The verification of the position of these parts of the curves is carried
out using three control test points defined in the European standards. [ H. Bromstrup, 2004]

In the  third regime  (C), at low  stress levels, the curve sharply  steepens, indicating  predominantly  brittle failure which  is
nearly load independent and is due to chemical aging and degradation of the polymer chains [B. Choi et al.,2009].

In simpler terms, regimes A and B depict "physical aging," while regime C signifies “chemical aging”. Regimes B and C
compete with each other. When chemical degradation occurs rapidly enough, regime B completely vanishes. Changes in the
molecular structure and morphology of materials, including factors like molecular mass, distribution, branch concentration,
and crystallinity, greatly influence crack initiation and PL. Enhanced polymerization processes and controlled adjustments
to these material parameters by suppliers have notably improved resistance to these forms of degradation [E.M. Hoang and
D. Lowe, 2008; P. Hutar et al., 2013]. Over the long service life of buried pipes spanning several decades, morphological
changes can occur, leading to physical aging. However, effective material stabilization can impede molecular degradation
and oxidation  processes. The type  and  concentration  of stabilizers  significantly  impact  the lifespan  of  pressurized pipes,
particularly in the quasi-brittle and brittle failure regions. Additionally, the concept of local crack tip aging helps explaining
the crack growth mechanisms associated with different stabilizer systems [G. Pluvinage and M. Elwany, 2007]

Figure.2. Schematic illustration of the failure behavior of water
pressurized PE pipe

Figure.3. Example of additional PL forces caused by the
existence of stone at the external surface of buried pipe.

Point Load

Traditional methods for estimating the residual lifetime of polymer pipes are mainly based on internal pressure tests defined
in  [EN ISO  9080, 2003] and  [ASTM  D2837, 2013]. However,  these  standard  procedures reveal  limitations,  particularly
concerning advanced materials like PE 100 or PE 100-RC which don't fail in these standard tests. Therefore, these methods
are not sufficient for predicting how long these materials will last.

To address these shortcomings, various accelerated tests like the NPT, PENT, FNCT, and CRB have emerged.[ N. Brown
and X. Lu,1991; M. Haager,2006] These tests aim at assessing the durability of polymer resins by evaluating their resistance
to slow cracking. The majority of articles propose a methodology to estimate the remaining lifetime of polymer pipes. They
utilize accelerated tests to measure creep crack propagation and apply fracture mechanics to describe cracks in the pipe wall
[P. Hutar et al., 2012; A. Frank et al., 2012]. However, the published articles predominantly concentrate on plastic pipes that
undergo only internal pressure testing [P. Hutar et al.,2012; A. Frank et al.,2012; L. Andena et al., 2009, E. Hoang and D.
Lowe, 2008; M. Farshad, 2004 ].

The increased SCG resistance of PE 100-RC presents difficulties for current short-term quality tests, as testing durations for
these materials  often surpass a year, leading to unacceptable delays and additional problems. For example,  in the PENT,

8
materials  no longer  develop brittle fractures  [A. Sukhadia,  2010] and the  detergent  Arkopal used  for the  Cone test  [ISO
13480, 1997] and FNCT degrades [F.L. Scholten et al., 2001].

In addition, there is criticism of the often-used FNCT due to the suboptimal sample shape (square instead of round), giving
rise to stress singularities. Furthermore, the multiple Round Robins show a large scatter between the different laboratories
[U. Niebergall et al., 2006], especially for long testing times.

Furthermore, as the usage of trenchless installation method for polymer pipe is increasing the influence  of surface damage
on the technical lifetime is becoming more and more important. So many testes have been designed for investigation of the
influence of PL on the lifespan of the polymer pipe.

There is therefore a need for new, quick and relevant tests. These have been found in the PL Test, used to determine the pipe
quality under point loading.

Point Load Test

Point loads are additional external loads concentrated on a small surface area. A schematic of additional PL forces caused
by  the  existence  of  stone  at  the  external  surface  of  buried  pipes,  is  shown  in  Figure.3.  The  PL  Test  was  developed  to
investigate the effect of this type of load, specifically to replicate the SCG failure caused by a rock pressing into an external
wall of polymer pipes, particularly in cases involving sand-less embedding of pipes for installation [PAS1075, 2009; 4 J.
Hessel, 2001]. The test is described in the standard [PAS 1075, 2009]. This test assesses the resistance of PE materials to
such damage, with PE100-RC being a high SCG-resistant PE developed for these applications. It is anticipated to exceed a
year in the PLT test. Additionally, an accelerated PLT+ test involving enhanced detergents and higher temperatures has been
introduced, with a minimum test duration of 450 hours for a PE100-RC rating [S. Nestelberger, J. Cheng, 2021].

Failure Mechanisms Under Point Loads

Buried pipes, including PE pipes, can experience undesired external forces, such as point loading, originating from rocks in
the soil or backfilling material. When a concentrated force is applied to the surface of the pipe, it creates localized stress on
the pipe wall  at the point of contact.  This load leads to permanent local stresses or deformations on the pipe wall.  These
additional stresses and deformations combine with the pressure-induced tensile stresses inside the pipe. The combined effect
of external point loading and internal pressure-induced tensile stresses can reduce the lifetime of buried pipelines, leading to
the initiation of cracks.

These cracks may start forming on the inner pipe wall, particularly near the area where the pipe is impacted by rocks. Once
initiated, it propagates from the inner side to the outer side of the pipe wall. The propagation of cracks is influenced by the
stress concentration at the point of loading and the material's response to these stresses. The propagation of cracks may result
in the formation of multiple macroscopic longitudinal lips along the pipe axis. These lips represent sections where the pipe
material has fractured, creating visible deformations along the length of the pipe.

The failure mechanism observed in these situations is described as brittle, indicating a lack of macro-ductility. This means
that the material does not undergo significant plastic deformation before failure, and the crack propagation is characterized
by  a  more  instantaneous  and  brittle  process.  This  suggests  that  the  failure  mode  is  similar  to  the  PL  process,  often
characterized  by  a  brittle  mode  of  failure.  This  brittle  mode  is  prevalent  in  low-pressure,  long-term  pipeline  service
operations. [ S. Nestelberger, J. Cheng, 2021; V. Rouyer and M. Cornette, 2001; J. lenz, 2001]

The studies highlight  the ongoing debate regarding the dominant mechanism  in buried pipes subjected to point loading is
stress relaxation or constant deformation.  Some  authors suggest  that in stiff  types of soil, the governing  process is stress
relaxation, indicating that pipes are more or less in a situation of constant deformation. On the contrary, other groups assume
a constant stress scenario, leading to a creep phenomenon. [V. Rouyer and M. Cornette, 2001]

The nature of the soil and the interaction with the buried pipe are still controversial topics. Different researchers have varying
perspectives  on  whether  stress  relaxation  or  constant  deformation  is  the  dominant  mechanism,  likely  influenced  by  the
specific properties of the backfilling materials and the surrounding soil. [V. Rouyer and M. Cornette, 2001]

Morphology of Cracks

To  improve  the  understanding  of  the  failure  process,  Figure.4.  illustrates  the  inner  wall  and  outer  wall  of  the  different
generation pipes during a PL test. When a dent presses into a PE pipe, it creates a dent on the outside and a bulge on the
inside of the pipe. As mentioned, the bulge on the inside leads to additional tensile stresses on the inner wall, in addition to
those caused by the internal pressure of the medium. The accumulation of tensile stresses, particularly at the bulge, can lead
to premature brittle fractures in the axial direction. The analyses of the crack in Figure 12-G shows that yielding occurs at

9
the center of the indentation, making the PE stronger due to molecular orientation. However, cracks initiate just next to the
center where no yielding occurs, forming a circle around the indentation. From the macroscopic cracking, the multiple cracks
around the indentation are described as brittle at the macro scale.

In laboratory testing PL, multiple cracking is observed due to very high-stress levels, whereas in practical applications, only
a  single  axial  crack  is  typically  observed.  This  difference  suggests  a  need  for  improvement  in  testing  methods  to  better
simulate real-world scenarios. [ E. van der Stok and F. Scholten, 2014]

A

D

G

B

E

H

C

F

I

Figure.4. Photograph of the HDPE pipes after PLT
A- B, C-D and E- F ) inner wall surface and the outer wall surface around the indentation for PE 80,
PE 100 and PE100-RC
G) for second generation PE pipe after PLT at 80℃ and 8 mm indentation in Dehyton show multiple cracking. The
crack starts just next to the center of the indentation where no yielding occurs[E. van der Stok and F. Scholten, 2014]
H and I) stamp on the outside of the PE 100 pipe and crack on the inner wall under the location of the
stamp [J. Lenz, 2001]

Effect of Various Parameters on The Point Load on The Lifetime of HDPE Pipes

The rocks can range in size and form, from small sharp-edged stones to large smooth  pebbles, to account for the diverse
range  of  rocks  encountered  in  trenchless  installations,  different  kinds  of  stamps  have  been  used  for  point  loading
measurements in the various researches.

Shape of Dentation

A few studies have examined the impact of dent geometries on the failure behaviors of polymer pipes. [Hutař et al., 2011].

[X. Zheng et al., 2019] present a novel approach using a notched HDPE ring specimen to assess the mechanical properties
of traditional notched HDPE pipes. Findings highlight reduced ultimate  loads with increased notch depth ratio, notably in
U-type and V-type grooves compared to L-type grooves. Additionally, the study suggests the potential substitution of NHP
specimens with  NHR ones for mechanical  property testing when the depth ratio is below 0.4 for different groove shapes.
The  proposed  finite  element  simulation  and  theoretical  model  successfully  predict  ultimate  loads,  aligning  well  with
experimental results obtained from tension tests.

Dimensions of Dentition

The modification of the stamp diameter exhibits that the stamp size has a crucial influence on the penetration depth [J. Lenz,
2001]. Study shows the penetration depth is smaller with a smaller stamp diameter than a larger one. While the penetration
depth did not even achieve 0.1 mm after 10 h with a small stamp diameter. The effects of the PL are higher with a bigger

10
stamp [J. Lenz, 2001]. the dent introduces additional tension stress on the inner pipe fibers and compression stress on the
outer area of the dent.  The degree of curvature in the dent directly influences the magnitude of these additional stresses. As
the dent becomes greater, the added stress intensifies. When you add radially acting compression stress, it forms an axial
compression stress state. Higher compressive loads mean higher reference stress and deeper stamp penetration. The curves
reflect this interaction between PL, dent, and resulting stresses [J. Lenz, 2001].

Outer Diameter

Different types of pipe failures, longitudinal, circumferential, or helicoidal, are determined mainly by the pipe diameter. In
small diameter pipes, where bending stresses are dominant, failures often  manifest circumferentially. Conversely, in large
diameter pipes, hoop stresses take precedence over bending stresses, resulting in longitudinal failure. If both bending and
hoop stresses are equally significant, the fracture path takes on a spiraled pattern [G. Pluvinage and M. Elwany, 2007]. This
phenomenon  is  observed  even  when  the  pipe  is  under  hydrostatic  pressure.  Consequently,  this  prompts  the  question  of
whether varying pipe sizes could influence stress distribution and loading conditions during the PLT.

with regard to the experimental results, it is evident that various pipe sizes exhibit unique deformation states when subjected
to a PL [S. Nestelberger, J. Cheng, 2021]. By maintaining the depth of the dent at a consistent 10% of the pipe wall thickness,
there is an escalation in stress intensity with an increase in pipe diameter. The higher stress intensity in larger diameter pipes
results in a quicker occurrence of SCG failure compared to smaller diameter pipes [Z. Zho and D. Chang, 2010]

Inner Pressure

The PL test was conducted with different hoop stress by varying water pressure between 1.5 MPa and 4.6 MPa. The results
illustrate a reduction in hoop stress corresponds to an increase in failure time. Reducing the internal pressure in the pipelines
in practice makes them less prone to failure due to point loading and helps extend the pipeline's lifetime We can also see that
decreasing the pressure by a factor of 2 (e.g., from 4 MPa to 2 MPa) extends the pipe's life by approximately 2.5 to 3 times.
[ E. van der Stok and F. Scholten, 2014]

2.3.

CHEMICAL DEGRADATION ASSESSMENT: EFFECTS OF CHLORINE

DISINFECTION

PE contains antioxidants, such as phenols, which can react with radicals. This suggests that disinfectants might compromise
the durability of PE pipes by destabilizing them, evident in the faster consumption of stabilizers, especially at the water-
polymer  interface.  This  reveals  how  disinfectants  could  impact  the  stability  of  the  pipes,  potentially  raising  durability
concerns. These chemical reactions alter the material's microstructure, which leads to property degradation and subsequently
causes the change of the failure mode.

These alterations, primarily attributed to oxidative degradation, include the formation of carbonyl groups, chain scission,
increased  crystallinity,  and  variations  in  melting  temperature.  Subsequently,  we  will  explore  how  these  alterations  can
influence mechanical properties.

Antioxidant Concentration Profiles

Antioxidants  play  an  essential  role  in  keeping  safe  polymers,  like  HDPE  pipes,  against  degradation  from  environmental
factors  such  as  oxidation.  However,  exposure  to  a  chlorine  medium,  like  chlorinated  water  or  water  containing  chlorine
dioxide,  significantly  reduces  the  concentration  of  antioxidants,  impacting  both  stability  and  mechanical  properties  [J.
Hassinen et al., 2004]. Chlorine dioxide seems to initiate an aggressive attack on antioxidants near the inner wall, employing
a single electron transfer process. This process extends deeper into the pipe, continually depleting antioxidants. The chemical
consumption surpasses 80% at the inner wall, indicating a substantial loss of stabilizer efficiency. [ J. Hassinen et al., 2004]

The  decline  in  antioxidant  concentration  makes  the  polymer  more  susceptible  to  degradation,  especially  from  oxidative
damage. This vulnerability alters the fracture behavior of the material. Weakened sections, particularly those closer to the
inner wall, are more prone to crack propagation. As a result, the mechanical strength of the material decreases, making it
more  brittle  and  less  resistant  to  stresses  and  pressures,  potentially  causing  ruptures  or  fractures.  This  heightened
susceptibility  may  lead  to  structural  failure  even  under  normal  operational  stresses.  When  comparing  the  effects  of
chlorinated water versus chlorine dioxide, the latter exhibits a notably faster rate of antioxidant consumption approximately
four times faster than chlorinated water [W. Yu et al. 2011]. This accelerated depletion further exacerbates the material's
vulnerability to oxidative damage, compromising its structural integrity and mechanical robustness.

Changes in Molecular Weight

11

Exposure of PE to chlorinated water can significantly impact its MW and distribution. Experimental results reveal that PE
oxidation is initiated in the presence of chlorine disinfectant. One primary aspect of this process involves the  buildup of
hydroperoxides which triggers chain scissions when their concentration is above a critical threshold. Reducing the length of
the polymer chains through chain scission results in a decline in the average MW and changes its distribution. This decrease
can lead to embrittlement, particularly when Mw gets close to approximately 70 kg /mol. Essentially, these changes in MW
can significantly affect properties such as thermal stability, ultimately impacting the performance and lifetime of PE pipes
in service. [ X. Colin et al., 2009] By Employing Saito relationship, the number of chain scissions 𝑠 and crosslink events 𝑥
per mass unit can be calculated from Mn and Mw variations.

1
𝑀𝑤
1
𝑀𝑛

−

−

1
𝑀𝑤0
1
𝑀𝑛0

=

𝑠
2

− 2𝑥

= 𝑠 − 𝑥

Eq.7.

Eq.8.

The competition between crosslink events and chain scissions during PE pipe exposure time is illustrated in Figure.5 [X.
Colin et al., 2007].

Figure.5. The number of chain scissions and crosslink events per mass unit calculated from the Saito equation against exposure time
by Colin et al

Crystallinity Profiles

The degradation of PE through chain scission leads to a decrease in its MW. The shorter chains can more easily incorporate
into the crystalline phase of the material, resulting in an increased crystallinity rate during aging. This phenomenon is known
as chemi-crystallization. At the nanometric scale, it is characterized by a reduction in the thickness of the amorphous phase
located between the crystalline lamellae. These microstructure modifications consequently contribute to the embrittlement
of  the  sample.  This  increase  in  crystallinity  was  measured  by  differential  scanning  calorimetry.  Figure.6.  illustrates  the
changes in the crystallinity ratio with the changes in the MW.

-
Figure.6. Crystallinity ration as a function of Mw
1/2 for PE

Fig.7.  Changes  in  the  ultimate  strain  𝜀𝑅  (■)  and  the
concentration of carbonyl groups (♦) during the exposure of PE
at 80℃

Influence on Mechanical Behavior

12

Chemical shifts, like increased carbonyl groups and chain scissions in consequence of changing in molecular mass, reduce
the material’s strength and ductility. Physical alterations, such as changes in crystallinity (change in morphology) and melting
temperature,  impact  flexibility  and  thermal  stability  [A.  Frank  et  al.2009]. Both  of  these  certainly  have  an  effect  on  the
mechanical properties of the pipes and on its resistance to crack initiation and SCG.

The  failure  mechanism  of  HDPE  moves  from  a  ductile  to  a  brittle  mode  as  the  corrosion  level  increases.  This  leads  to
subcritical crack propagation, which deteriorates the load capacity of the structure. Field observations and pressure testing
have demonstrated that exposure to chlorinated water leads to premature brittle fractures in PE pipes. [J. Dear and N. Mason,
2001; J. Sanders et al., 2009]

[B. Fayolle et al., 2007]. demonstrated that, in PE strain at break decreases with carbonyl buildup, as shown in  Figure.7.,
illustrating strain at break against exposure time at 80 °C, including the kinetic curve of carbonyl buildup. The observation
is  that  embrittlement,  signified  by  a  catastrophic  decrease  in  strain  at  break,  occurs  when  the  carbonyl  concentration
approaches 0.1 mol kg-1.

Upon comparing Figure 17 and 18, it is evident that carbonyl concentration aligns with the growth in the number of chain
scissions and crosslink events, suggesting a correlation between the magnitude of carbonyl concentration and the number of
chain scissions, as observed in earlier studies [B. Fayolle et al, 2007]

Fractography

Figure 9 [W. Yu et al., 2011] shows crack that was formed during hydrostatic pressure testing of a pipe exposed to water
containing chlorine dioxide. Around the major crack, there are areas of brittle material indicated by the presence of small
cracks. These probably indicate parts of the material that are highly degraded. Interestingly, the primary crack advanced into
a region of fresh material but ceased further propagation. Then, some aggressive chemicals attacked the material near the
crack, suggesting the potential for continued crack growth.

This mechanism, we can name "degradation-assisted crack propagation", explains the possible causes for premature fractures
in even thick-walled pipes, even if only the surface of the material is affected by degradation. Figure10 illustrates the different
phases involved in degradation-assisted crack growth, providing a visual representation of this phenomenon [W. Yu et al.,
2011].

Figure.9. Scanning electron micrograph of a crack in a pipe
after 121h of exposure to water containing 4ppm chlorine
dioxide.

Figure .10. Schematic representation of degradation assists
crack propagation. The arrows indicate the principal stress
direction

2.4.  LIFETIME PREDICTION: COMBINED INFLUENCE OF PUNCHING AND

CHEMICAL DEGRADATION

The  study  of  the  combined  influence  of  mechanical  stress  and  chemical  degradation  on  the  lifetime  of  polymer  pipes
represents a critical unexplored area within  the current scientific literature. Particularly, there is a notable absence of prior
research addressing the combined effect of chlorinated water exposure and point loading on the assessment of polymer pipe
durability, signifying a substantial gap in existing knowledge.

In terms of empirical findings, there hasn't been much information about how mechanical and chemical factors work together
to influence the lifetime of polymer pipes. However, [A. Tripathi et al., 2021] took a different approach by using a computer
model to predict this combined impact. They employed a coupled chemo-mechanical modeling approach to simulate stress

13
corrosion cracking in HDPE exposed to a bleach solution. The coupled chemo-mechanical model is implemented in a finite
element simulation using ABAQUS.

The model introduces a simplified corrosion model, specifically focusing on tracking the diffusion and reaction of HOCl.
Application of the model to predict the behavior of HDPE specimens exposed to a corrosive environment under continuous
loading  produces  stress-life  curves,  As  shown  in  Figure.11.  The  simulations  reveal  two  distinct  regimes  for  unexposed
specimens and three regimes for exposed ones, demonstrating the response of material to different stress levels and corrosion
effects.  Notably,  in  the  high-stress  regime,  corrosion  has  a  negligible  impact  on  failure,  with  yielding  dominating.
Conversely, in the low-stress regime, corrosion becomes a significant factor, leading to failure through craze breakdown.

Figure.11. Simulated stress-life curves.

In  this  model,  the  fracture  kinetics  of  HDPE  have  been  predicted  under  different  bleach  concentrations,  offering  a
comprehensive  understanding  of  crack  growth  rates.  The  dependency  of  crack  growth  rate  on  bleach  concentration  is
characterized  by  two  regimes:  a  constant  rate  influenced  by  corrosion  and  a  power-law  model  representing  mechanical
loading. The outcomes of this research show that this predictive model emerges as a valuable tool for estimating the service
lifetime  of  HDPE  structures.  It  uniquely  considers  the  interplay  between  chemical  exposure  and  mechanical  stress,
addressing the essential aspects overlooked by existing studies.

It is important to note that they simulated a simple tested square rather than a pipe. This approach could limit the reliability
of the results in the real world. Additionally, they used HOCl as the chemical species responsible for the degradation, but
this cannot be the case since it is not a radical. Instead, it is HO° and Cl° that attack the antioxidants and PE. This constitutes
another limitation to consider. The exploration of the combined influence of punching and chemical degradation on polymer
pipe lifetime  prediction is yet unknown,  and it holds practical implications  for real-world  performance. By doing so, the
project seeks to provide a comprehensive understanding of the synergistic effects that these factors may exert on the material
properties and structural integrity of polymer pipes.

3. Scientific Approach

3.1.  MAIN IDEAS

The central objective of this project is to comprehensively understand the premature failure mechanisms of PE pipes under
point  loading,  exploring  the  interplay  between  mechanical  stresses  and  chemical  degradation  induced  by  chlorinated
disinfectants. This investigation  aims to provide insights  into how these factors jointly influence pipe failure and develop
predictive models for estimating pipe Lifetime.

To achieve this goal, the scientific approach is divided into three phases: 1. examining the impact of punching on mechanical
stress and strain distribution within the pipe via the development of finite element software, 2. assessing chemical degradation
induced by chlorine disinfection and 3. predicting the pipes'  lifetime by considering the combined influences of punching
and chemical degradation.

Phase 1: Computational Modelling

In this research phase, we will use the finite element software to create precise models of PE pipes and simulate various PL
scenarios accurately. By adjusting factors like the shape and dimensions of the punch, the punching depth, the dimensions
of the pipes, and the internal water pressure, we can replicate real-life mechanical pressures on the pipes.

The FEA analysis will compute how stress and strain are distributed within the pipes' structure, particularly focusing on the
initiation and propagation of fractures caused by point loading. These findings will offer crucial insights into how well PE

14
pipes withstand mechanical pressure and will aid in estimating their durability under different stress conditions. Additionally,
our parametric analysis  enables us to identify  critical damage  conditions, exploring various parameters to understand the
thresholds  for  damage  occurrence.  Furthermore,  we  aim  to  propose  effective  solutions  aimed  at  reducing  or  preventing
damage  through  innovative  approaches. This involves  damage  modeling,  seeking  to develop comprehensive  models  that
depict the behavior of the pipes under diverse stress situation.

Phase 2: Chemical Degradation Assessment: Effects of Chlorine Disinfection

In the second phase, the project will delve deeper into the effects of chlorinated disinfectants on the stability and longevity
of PE pipes. The primary objective is to predict how various levels of chlorination affect PE pipes over time. This will be
achieved by observing changes  in material  properties and degradation mechanisms.  This phase involves  conducting  both
physical and chemical characterization tests on materials, as well as assessing mechanical properties at different levels of
degradation under chlorine exposure

Phase 3: Lifetime Prediction: Combined Influence of Punching and Chemical Degradation

In  the  final  phase  of  our  research,  we  will  combine  the  findings  from  the  study  of  mechanical  stresses  and  chemical
degradation to predict the lifetime  of PE pipes more  accurately.  The aim  is to understand how  mechanical  and chemical
factors together affect the integrity  of the pipes. We will  use the data from  both these analyses  to create models that can
simulate the combined effects of point loading and chemical exposure over time. We will also study how continuous chemical
attack affects the material properties of the PE pipes when they are subjected to mechanical point loading.

A  multi-Physical  Problem  Solving  by  solving  for  the  interaction  between  chemical  and  mechanical  degradation  will  be
carried out to identify which factors most significantly reduce the pipes' service life and under what conditions. The predictive
model developed will aim to forecast the point of failure more accurately, allowing for early interventions to improve  pipe
durability and reliability.

3.2.  THE TOOLS AND THE TECHNIQUES

In this section, we will describe the analytical equipment and protocols used for each phase. I will detail also the technical
tools employed.
Phase 1: Computational Modelling

For modeling our scenario, Abaqus/CAE 2022 teaching version was used to examine the stress/strain situation in the pipe
wall during PLT. The objective is to create accurate models that can investigate the stress/strain distribution in the pipe wall
throughout  each  time  step  of  the  loading  sequence  (punching  –  heating  –  internal  pressure).  Additionally,  it  enables  the
calculation of results such as the required force during indentation or the resulting wall thickness reduction below the pin.

The methodology is divided into two parts with different steps. Table 1 summarizes the scientific method used in this phase,
each focusing on different aspects of the modeling process.

Part 1: Initial Modeling and Validation

Pipes under internal pressure: First, I modeled a general case of pipes under internal pressure (diameter extern 32mm with
thickness  2.9mm  under  10bar  internal  pression),  creating  both  3D  and  2D  models.  To  optimize  computational  time  and
resources, the models were simplified. Six different cases were created: full pipes, half-pipes, and quarter-pipes, each in both
3D  and  2D.  The  ABAQUS  results  were  then  compared  with  corresponding  analytical  results  (eq.  1),  focusing  on
circumferential stress, as these stresses are more critical than radial or longitudinal stresses. Once this step was validated,
another approach was taken to further validate our model.

Pipe under Point Load: In this step, we began by selecting a reference article to guide our modeling efforts. The scenario
described in the reference article was modeled to identify optimal parameters for our simulations in Abaqus, such as finding
the best options for contact between the pin and pipe. In the article, a pipe diameter of 32 mm, SDR 11, with a punch diameter
of 10 mm and depth of punching equals 8% of diameter extern of pipe were simulated. This stage provided a baseline for
parameter identification in the simulation with Abaqus. For the material properties, we used the values from the article, but
many parameters were not mentioned. Therefore, we searched other articles, handbooks, or data sheets of PE100 pipes to
find values that are adequate to reality.

Adding Thermal and Pressure Loads: In the third step, after validating the model of the pipe under point load, thermal and
internal pressure loads were applied to the models step by step. Specifically, a temperature of 80°C and an internal pressure
of 885 bar (equivalent to 4 MPa hoop stress) were applied.

15
Elasto-Plastic  Modeling:  In  the  absence  of  specific  details  on  elasto-plastic  studies  from  the  reference  article,  our  study
focused on exploring various plasticity models available in Abaqus. Understanding how materials transition from elastic to
plastic  behavior  is  crucial,  and  we  utilized  the  Von  Mises  criterion  to  define  the  yield  surface,  a  fundamental  aspect  in
material modeling.

The Von Mises yield surface is a cylinder in principal stress space, defined by the following formula:

𝑓 ( 𝜎𝑖𝑗 , 𝜀̅𝑝 ) =   √

1
2

[(𝜎1 −   𝜎2)2 + (𝜎1 −   𝜎3)2 + (𝜎2 −   𝜎3)2  ]    −  𝑌(𝜀̅𝑝) = 0

Eq.9

The radius of this cylinder varies with the yield strength, with its axis parallel to the line σ1=σ2=σ3. If the stress state  is
inside the cylinder (f < 0), the material behaves elastically. If the stress state is on the cylinder surface (f = 0), the material
begins to deform plastically. The yield strength 𝑌(𝜀̅𝑝) can increase during plastic deformation, indicating that 𝑌  is a function
of the total plastic strain 𝜀̅𝑝.

Abaqus provides several hardening models: isotropic, kinematic, and combined hardening. Each model influences how
materials respond under varying loading conditions. Figure.12. illustrates the various hardening models in Abaqus.

•

Isotropic Hardening: This model assumes uniform hardening in all directions, leading to equal resistance to
plastic deformation throughout the material. It simplifies loading and unloading scenarios by maintaining
consistency in the material's response.

•  Kinematic Hardening: More intricate, this model accounts for the movement of the yield surface. It effectively
models phenomena such as the Bauschinger effect, where material behavior differs between loading and
unloading cycles. This model is particularly suitable for simulations involving cyclic loading conditions.

•  Combined Hardening: This model combines aspects of isotropic and kinematic hardening, offering a more

comprehensive representation of material behavior under complex loading scenarios. It proves beneficial for
materials subjected to varied and cyclic loading conditions.

Figure.12. The various hardening models in Abaqus: Isotropic Hardening(left), Kinematic Hardening (middle), Combined
Hardening (Right)

Since our model does not involve cyclic loads, we opted for the isotropic hardening model in our simulations. However,
we acknowledge that the yield strength may vary with plastic strain. In isotropic hardening, yield strength can be in three
forms: perfect plasticity (constant 𝜎𝑌), linear strain hardening (𝜎𝑌 increases linearly with 𝜀̅𝑝), and power law hardening
(where 𝜎𝑌 follows a power law relationship with 𝜀̅𝑝).

16

𝜎𝑌 = 𝑐𝑜𝑛𝑠𝑡𝑎𝑛𝑡

𝜎𝑌 (𝜀̅𝑝) =   𝑌0 + ℎ 𝜀̅𝑝

𝜎𝑌 (𝜀̅𝑝) =   𝑌0 + ℎ (𝜀̅𝑝)1/𝑚

Figure 13. the different form of yield stress

By evaluating these models, we aimed to understand their impact on simulation results and their relevance to real-world
material behavior. Our approach underscores the importance of selecting an appropriate hardening model to achieve
accurate simulations, ensuring that our results align closely with practical observations and theoretical predictions in
material science and engineering.

Step  Description

Parameters

objective

General modeling of pipes under
internal pressure and comparison of
ABAQUS results with analytical
results
Modeling specific scenarios based on
reference article
Application of heat and internal
pressure

3D and 2D models of full, half, and
quarter pipes; comparison of
circumferential stress with analytical
results
Diameter: 32 mm, SDR 11, Punch: 10
mm diameter, 8% depth for PE100
Heat: 80°C, Pressure: 885 bar (4 MPa
hoop stress)

Elastic and elasto-plastic modeling

Different modeled hardening

1

2

3

4

Optimization of computational resources
and validation of model accuracy

Data retrieval and parameter uncertainty

Verification against reference results

Matching model parameters with
reference data

Table.1. Summery of the method scientific in the first phase

Part 2: reel case and parametric study

Once our model is validated through comparison with reference data, we can proceed to the real case study. We will utilize
the same modeling parameters but with different material properties and boundary conditions. Specifically, we will model
¼ of the PE100 pipe with an outer diameter of 90 mm and a wall thickness of 7.8 mm. The pipe will be subjected to a punch
with a 10 mm diameter, penetrating to a depth of 8.2% of the pipe's external diameter. The system will be under a pressure
of 6 bar at a temperature of 40°C.

The mechanical properties of the PE100 material used in this study are derived from the tensile tests we have conducted.
These properties include Young's modulus, yield strength, tensile strength, and elongation at break, which are essential for
accurately simulating the pipe's behavior under mechanical stress.

A key objective of this study is to identify critical points where failure is most likely to occur. This involves determining the
stress  concentrations  and  analyzing  the  failure  modes  of  the  pipes  under  purely  mechanical  loading  conditions.  By
understanding the stress thresholds at which failure becomes imminent, we can gain insights into the material's performance
and potential failure mechanism.
Phase 2: Chemical Degradation Assessment: Effects of Chlorine Disinfection

In  this  phase  of  the  project,  we  focus  on  characterizing  the  materials  both  before  and  after  degradation  due  to  chlorine
exposure. This ensures a comprehensive understanding of material degradation and performance under real-world conditions.
To  achieve  this,  a  series  of  tests  were  initiated  to  characterize  the  physical  and  chemical  properties  of  the  materials.
Additionally,  mechanical  tests  were  conducted  to  assess  the  material  properties  under  stress.  To  analyze  changes  in  the
properties of materials at different stages of degradation, we prepared a platform to expose the materials to a chlorine solution

17
(with a concentration of 4 ppm HoCl at 40°C). Material samples are collected monthly to analyze changes in their properties
at different stages of degradation. This phase is divided into three main parts: Film Fabrication Using Compression Molding,
Physico-Chemical Characterization of Materials, Mechanical Behavior of Materials.

The material studied

The material used in this study is high density polyethylene, HDPE (PE100) used in the drinking water pipes. It was then
extruded to make stoppers. It is a semi-crystalline thermoplastic, comprising an amorphous phase and a crystalline phase, in
the form of a spheroidic aggregate. The PE100 grade is one of the latest developments in polyethylene resins, and as a result,
guarantees the best mechanical characteristics of all polyethylene resins on the market.

Part 1: Film Fabrication Using Compression Molding

To ensure uniform degradation throughout the material and avoid any gradients, specimens with a thickness of less than 1
mm were prepared. This is crucial for accurately assessing the material properties at various levels of degradation, it was
decided to produce a film with a thickness of 600 microns. Consequently, this process was divided into three steps:

Step 1: Characterize materials before film fabrication: Initially, we performed physical and chemical characterization of the
materials before starting to fabricate the  films under pressure  and compression.  This  was done  to determine  whether the
properties would change after the manufacturing process. For this reason, we have done the ATG, FTIR, OIT, and DSC. we
did also the mechanical tensile test the pipes.

Step 2: Fabricate films under various conditions using a compression molding machine: We worked on PE100 pipes so
These  pipes  were  cut  and  ground  into  small  granules  suitable  for  film  fabrication  using  a  Gibitre  compression  molding
machine. Our goal was to produce thin films free of defects, with a smooth and flat surface, by varying the conditions of
temperature, pressure, and molding time. The following table summarizes the conditions tested:

Film

1
2
3
4
5
6
7

Temperature

Molding Time

Pressure Time

Cooling Condition

(°C)
160
160
160
160
160
180
180

(s)
120
120
150
120
120
120
150

(s)
120
60
60
90
120
120
120

Fast
Fast
Fast
Fast
Slow
Slow
Slow

Table.2. Summary of the conditions used in fabrication the films

 Step 3: Re-characterize films to assess changes in material properties due to fabrication conditions: After fabricating the
films, we repeated the physico-chemical tests (TGA, FTIR, OIT, DSC) to assess any changes in material properties due to
the  fabrication  process.  This  re-characterization  was  essential  to  determine  if  the  material  properties  had  altered  under
different fabrication conditions.

Finally, we decided to fabricate the films at the temperature of 180°C and waiting 150s for molding the material et after that
press the materials under pressure 22 bar for 120s because the films produced under them had a very flat shape without any
signification changes in the properties. in additional we have chosen the slow cooling condition because Rapid cooling was
found to decrease the degree of crystallinity.

Part 2: Physical-Chemical and Mechanical Behavior characterization before/ after degradation

This  section  presents  the  study  of  the  physical-chemical  and  mechanical  characterization  of  materials  before  and  after
degradation. The study is divided into two main steps:

Step  1:  Initial  Characterization:  After  fabricating  the  films,  materials  without  degradation  were  characterized  using  the
aforementioned tests to establish a baseline.

Step 2: Re-characterization: All characterization tests were then conducted on samples that aged over time. Material samples
were collected monthly to analyze changes in their properties at different stages of degradation.

18
For this purpose, various tests  were  conducted, including  TGA (Thermogravimetric Analysis), FTIR (Fourier Transform
Infrared Spectroscopy), OIT (Oxidation Induction Time), DSC (Differential Scanning Calorimetry), and tensile tests. In the
following sections, I will explain in detail the methodologies and equipment used for these tests.

Mechanical characterization: uniaxial tensile tests: In order to assess the mechanical properties of the materials, uniaxial
tensile tests were systematically conducted. Initially, mechanical tests were performed on samples procured from Suez to
establish their baseline mechanical characteristics. The samples underwent precision polishing to achieve the requisite shape
for testing. Subsequently, thin specimens were fabricated with a thickness of approximately 6 microns, specifically designed
for mechanical evaluation.

Figure 14. Specimens for Tensile Test in Two Thicknesses and the Testing Machine

Following fabrication, tensile tests were performed on these thin specimens’ post-degradation to investigate how exposure
to chlorine and other degradative factors influenced their mechanical behavior over time. This step aimed to provide insights
into the material's durability and performance under simulated environmental stressors.

Given  the  operational  context  where  pipes  are  exposed  to  a  40°C  environment,  the  tensile  tests  were  conducted  at  this
temperature to simulate realistic operating conditions. During testing, specimens were subjected to controlled displacement
until failure, with meticulous recording of stress-strain data.

The testing procedures adhered to relevant standards and norms to ensure consistency and reliability of results. Tensile tests
were conducted using an Instron 5966 machine equipped with a precision 10 kN load cell. To account for variability, ten
tests were performed on samples in their initial state, with an equivalent number conducted on aged samples. The test speeds
were meticulously set at 10 mm/min for thicker samples and 5 mm/min for thin films, focusing on assessing both breaking
behavior and elastic properties in accordance with ISO 527 standard guidelines.

Thermogravimetric Analysis (TGA): To determine the thermal stability and composition of the materials, we measure
weight changes as a function of temperature, including the mass of carbon black. Samples were heated at a controlled rate
in an inert atmosphere  with two different  gas, and  weight  loss  was recorded. In this study  TGA of studied samples  was
carried out using a TA Instruments Q500 analyzer. At first cycle the samples were heated from temperature 30°C to 600°C
at a heating rate of 20°C/min under Nitrogen et in the second cycle they were heated from temperature 600°c to 700°C at a
heating rate of 20°C/min under oxygen.

FTIR Spectroscopy. To observe changes in a material caused by disinfectant chlorine, HDPE/PE100 was characterized by
Fourier transform infrared spectroscopy (ATR FTIR) in the range from 4000 to 650 cm−1. Each sample was pressed with the
flat  pressure  tip  at  the  maximum  pressure.  From  each  sample,  8  parallel  spectral  measurements  were  carried  out.  We
Identified functional groups and molecular structures present in the materials before and after degradation.

Oxidation Induction Time (OIT): We effected the OIT to measure the oxidative stability of the materials. Samples were
heated in a  differential scanning calorimeter (DSC)  on  nitrogen  rich environment  with  flow rate  50  ml/min  at  first pour
arrived a temperature stable et after that we change le gas to the oxygen for a period more than 3h for recording the time to
oxidation onset. We conducted at two different temperatures, 220°C and 200°C, with the following conditions for each:

Méthode1:

➢  Selected Gas N2
➢  Ramp 10°C/min to 220°C
➢
Isotherm 5min
➢  Selected Gas O2
➢

Isotherm 360min /240min

Method 2:

➢  Selected Gas N2
➢  Ramp 20°C/min to 200°C
➢
Isotherm 5min
➢  Selected Gas O2
➢

Isotherm 360min /240min

Table 3. The methodologies used for OIT test

19

Differential Scanning Calorimetry (DSC): DSC studies were performed on a DSC 25 - TA Instruments (TA Instruments).
Samples between 5.0±0.01 et 10.0±0.01 mg were sealed in aluminum pans, heated from 25°C to 220°C at heating rate of
10°C/min. Nitrogen flow rate was 40 ml/min. Measurement of each sample was performed three times. All experiments were
carried out under nitrogen atmosphere.

We determined the thermal transitions such as melting and the rate of crystallinity. Samples were heated at a controlled rate,
and heat flow was measured to identify thermal events. Rate of crystallinity (χc) was calculated using the formula. where
∆H°f is the enthalpy of fusion for 100% crystalline polyethylene (290 J/g).

𝜒𝑐  =   ∆𝐻𝑓 / ∆𝐻°𝑓,

Eq 10.

Phase 3: Lifetime Prediction: Combined Influence of Punching and Chemical Degradation
The third phase of the project aims to predict the lifetime of polyethylene (PE) pipes when subjected to both mechanical
punching  stresses  and  chemical  degradation.  This  involves  a  detailed  analysis  using  experimental  data  and  advanced
simulation techniques.

The primary objective of this phase is to quantify how chemical degradation, induced by chlorinated substances, exacerbates
the failure of PE pipes under  mechanical loading. By understanding the  compounded effects of these  factors,  we aim to
elucidate their role in the premature failure of PE pipes.

This advanced approach will simulate the internal layer degradation and its impact on the overall pipe structure when exposed
to both mechanical stresses and chemical attacks from chlorinated disinfectant

Methodology

Although Phase 3 has not yet been conducted, the planned approach involves the following steps:

At first, data collected in Phase 2, which includes the physico-chemical characterization and mechanical properties of the
degraded materials, will be utilized. This data will help in creating accurate models of degraded PE pipes.

After that, the internal layers of the simulated PE pipe models will be modified to reflect a reduction in mechanical
properties due to chemical attack, based on the degradation data obtained from Phase 2. Using ABAQUS, we will
dynamically alter the properties of the PE pipe to model real-time degradation effects occurring under service conditions.
This involves:

▪  Establishing Baseline Mechanical Properties: Initial mechanical properties of the PE pipes will be determined.

▪

Introducing Chemical Degradation Parameters: Chemical degradation parameters will be introduced to simulate
the weakening of the material over time.

▪  Applying Mechanical Loads: Mechanical loads will be applied to observe the interplay between punching

stresses and the chemically weakened state of the pipes.

4. Results and discussion

After describing the study material, test protocols, and analysis techniques used, we will now focus on the results of our
study, starting with the firs phase of the project.

4.1.  PHASE 1: COMPUTATIONAL MODELLING

These simulations allow for precise analysis of the  mechanical behavior of the material under various senario. However,
accurate  results  are  contingent  upon  correct  model  calibration.  Therefore,  in  the  initial  part  of  the  modelling  phase,  we
verified and validated our model with different approaches. This validation process was divided into four sections, and the
following are the results for each part:

Part 1: Initial Modeling and Validation

Pipes under internal pressure: Here, we present the results of the first part of our study, which examines a pipe with an
external diameter of 32 mm and a wall thickness of 2.9 mm under an internal pressure of 10 bar. We compare the results
from six models of pipes in 2D and 3D configurations (complete pipe, half pipe, and quarter pipe) using Abaqus simulations
and analytical calculations based on Equation 5. In the figure.15, we have plotted the circumferential stress distribution from
various simulation methods as a function of the distance between the inner and outer radius, as this is the most significant
stress in the case of internal pressure.

20

)
a
P
M

(

n
o
i
t
u
b
i
r
t
s
i
d

s
s
e
r
t
s

l
a
i
t
n
e
r
e
f
m
u
c
r
i
c

5

4.5

4

3.5

3

2.5

2

1.5

1

0.5

0

Analytique
Pipe_Quart_3D
Pipe_Quart_2D

Pipe_Complet_3D
Pipe_Complet_2D

Pipe_Half_3D
Pipe_Half_2D

13.1

13.6

14.1
14.6
15.1
The Path between Ri and Re (mm)

15.6

16.1

Figure 15. the circumferential stress distribution from various simulation methods

As shown, the results from the 2D simulation exactly match the analytical results. However, there are noticeable differences
when comparing the 3D simulation results to both the 2D simulation and the analytical results. This finding underscores the
utility of 2D simulations for achieving results that are exactly matched to theoretical predictions that in the specific case of
a pipe under internal pressure, simplifying the model from 3D to 2D does not compromise accuracy. The circumferential
stress,  being  the  most  critical  stress  under  internal  pressure,  is  well  captured  by  the  2D  simulations.  This  simplification
significantly reduces computational time and resources while maintaining the fidelity of the results

the discrepancy can be attributed to the limitations of the student version of the software, which restricts the mesh refinement.
Consequently, the 3D model exhibits larger deviations. Nevertheless, it is evident that with a highly refined mesh, the 3D
simulation  results  can  closely  align  with  the  analytical  predictions.  Additionally,  by  optimizing  the  mesh  size  in  3D
simulations, we can achieve results that are nearly identical to those obtained analytically. Therefore, the choice of mesh size
is crucial for ensuring accurate simulation outcomes. This study demonstrates that, under the conditions tested, a 3D case
can effectively be simplified to a 2D case without altering the results. This simplification is particularly beneficial in practical
engineering applications where computational efficiency and resource optimization are paramount.

21

5

4.5

4

3.5

3

2.5

2

1.5

1

0.5

0

)
a
P
M

(

n
o
i
t
u
b
i
r
t
s
i
d

s
s
e
r
t
s

l
a
i
t
n
e
r
e
f
m
u
c
r
i
c

Analytique

Pipe_Complet_3D

Pipe_Complet_3D -mesh Fin

13.1

13.6

14.1
14.6
15.1
The Path between Ri and Re (mm)

15.6

16.1

Figure 16. The effect of the mesh size in the resultant.

Pipe under Point Load: In the second part of our study, we modeled a pipe with an outer diameter of 32 mm and a thickness
of 2.9 mm subjected to a point load applied by a 10 mm diameter pin, penetrating 8% of the pipe's outer diameter at the
contact point. This modeling was conducted in 2D using PE100 material properties: a Young's modulus of 950 MPa and a
Poisson's ratio of 0.35, assuming elastic behavior.

Initially, we attempted to simplify the problem by modeling it in 2D. However, we encountered convergence issues when
comparing the results with reference values. By switching to a 3D model with the same parameters, we achieved results that
aligned  more  closely  with  the  reference  data.  This  indicated  that  for  point  load  cases,  a  2D  model  was  inadequate.
Consequently, we proceeded with 3D modeling for subsequent analyses.

In the 2D model, we experimented with different parameters to approach the reference results, such as the position and shape
of the plateau, refining the mesh, and varying the contact parameters between the pin and the pipe, including the coefficient
of friction. Despite these adjustments, significant improvements were not observed. Conversely, when modeling the same
problem in 3D with identical parameters, the results varied substantially

The circumferential strain distribution under selected loading conditions, as shown in Figure 18., was taken from the inner
fiber from the pin to the end cap (indicated by the red arrow in Figure 17).

Figure 17. Simulated of a 32x2.9mm pipe during PLT condition (8% deflection of OD); arrow indicates inner fibre position.

0.084587554

0.07683

0.0502

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0.09

0.08

0.07

0.06

0.05

0.04

0.03

0.02

0.01

0

22

Article Referance

Abaqus -M25000

Abaqus -MF-1000

Abaqus -M300

Abaqus-M90

Abaqus-2D

0

5

10

15

20

25

30

35

Inner fiber from the pin to the end [mm]

Figure 18. Circumferential strain computed along the inner fiber of the pipe wall starting below pin position

The  inadequacy  of  the  2D  model  for  the  punching  case  can  be  attributed  to  its  inability  to  capture  the  complex  three-
dimensional stress distribution around the punch. Significant out-of-plane stresses and deformations arise during punching,
which a 2D model cannot accurately represent. In contrast, the 3D model effectively captures these complexities, results in
better agreement with the reference results.

In the 3D case, we observed that refining the mesh improved the results towards the reference values.  However, we also
noted that beyond a certain level of mesh refinement, the results did not change significantly, indicating mesh convergence.
This  effect  is  illustrated  in  the  figure  19.  Therefore,  there  exists  an  optimal  mesh  size  for  the  model,  which  we  refined
primarily around the point of punching.

0.09

0.08

0.07

0.06

0.05

0.04

0.03

0.02

0.01

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0

0

5000

10000

15000

20000

25000

Number of meshes

Figure 19.the effect of convergence the mesh

Additionally, the observation of mesh convergence underscores an important aspect of finite element analysis. While refining
the mesh can enhance accuracy, there is a point of diminishing returns where further refinement does not significantly alter
the  results.  Understanding  this  convergence  behavior  is  crucial  for  conducting  efficient  and  reliable  finite  element
simulations.

23
We  also  found  that  by  adjusting  certain  parameters,  we  could  improve  our  results.  For  instance,  using  the  "Not  Allow
Separation" option in combination with "Hard Contact Normal Behavior" ensured that the contact surfaces, once engaged,
did not separate during the analysis. This setting is particularly important in simulations where a permanent bond is needed,
such as with a pin inserted into a pipe where separation must be avoided once contact is established.

Article

Abaqus -M25000-Not allow Sep

Abaqus -M25000

0.084587554

0.0793055

0.0827315

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0.09

0.08

0.07

0.06

0.05

0.04

0.03

0.02

0.01

0

0

5

10

15

20

25

30

35

Inner fiber from the pin to the end [mm]

Figure 20.the effect of Not Allow Separation in the results

In table 4. you can see the parameter that we chose for the following steps:

Category
Geometry

Material

Interaction

Load

Mesh

Parameter
Pipe
Pin Size
Support
Type
Elastic Modulus (E)
Poisson's Ratio (ϑ)
Contact Type
Normal Behavior
Tangential Behavior
Depth of Indentation
Displacement for the Dent (U_y)
Mesh Size Distribution

Value
32 x 2.9 mm (SDR11)
10 mm (Standard 1075 PLT)
60 mm with angle 120°
PE100
950 MPa (20°C)
0.35
Surface to surface contact
Hard contact, no separation
Penalty, friction 0.3
8% of D_ext
-2.56 mm
Biased, high accuracy near the dent

Table1. Parameters used in the 3D PLT model

Adding  Thermal  and  Pressure  Loads:  According  to  the  referenced  article,  the  Young's  modulus  (E)  for  the  material  in
question  changes  from  950  MPa  at  20°C  to  300  MPa  at  80°C.  Despite  this  adjustment,  our  observations  indicated  that
circumferential deformation did not significantly change, as illustrated in Figure 21.

The initial observation that circumferential deformation remained unchanged suggests the need to consider other deformation
modes  to  accurately  capture  the  material  behavior  under  thermal  loads.  Recognizing  that  mechanical  deformation
encompasses both mechanical and thermal components, we included thermal deformation in our analysis. This inclusion
proved crucial, as it resulted in a deformation change of at least 20 %, demonstrating that thermal deformation is a significant
factor. The referenced article and data sheet indicate that the coefficient of thermal expansion (α) for PE100 ranges between
1.3×10⁻⁴ K⁻¹ and 5.1×10⁻⁴ K⁻¹, varying with temperature.

0.12

0.1

0.08

0.06

0.04

0.02

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0

0

5

24

Article Referance
Modelisation step ppoint load
Modelisation point load + chaleur

10
25
Inner fiber from the pin to the end [mm]

15

20

30

35

Figure21. the effect of changing young modulus in the Circumferential strain computed along the inner fiber of the pipe wall starting
below pin position after heating until 80°C-

To address the impact of the coefficient of thermal expansion and to determine the appropriate value used in the referenced
article for the subsequent validation steps, we conducted a parametric study. Ultimately, we found that a value of 1.5 yielded
nearly perfect results. However, from the results, we can conclude that for a more realistic scenario, it is preferable to account
for the temperature dependence of the coefficient of thermal expansion. Figure 22 demonstrates the impact of the coefficient
of thermal expansion as reported in various studies.

The variability of the coefficient of thermal expansion with temperature further complicates the analysis. As temperature
increases, the coefficient of thermal expansion increases, leading to an augmentation in the overall thermal deformation. This
nuanced understanding is crucial because it underscores the  importance of incorporating temperature-dependent  material
properties into finite element models to achieve accurate predictions.

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0.16

0.14

0.12

0.1

0.08

0.06

0.04

0.02

0

0

Article Referance
𝜶 (𝑻) : 2,2 to 3,5 *10-4
𝜶 -2*10-4
𝜶 : 2,2*10-4
𝜶 : 3,5*10-4
𝜶- 1,5*10-4

25

30

35

5
Inner fiber from the pin to the end [mm]

10

15

20

Figure22. The effect of the thermal expansion coefficient in the Circumferential strain computed along the inner fiber of the pipe wall
starting below pin position after heating until 80°C

From  the  preceding  steps,  we  selected  a  coefficient  of  thermal  expansion  of  1.5×10⁻⁴  K⁻¹.  We  then  added  0.885  bar  of
pressure  on  the  internal  walls  of  the  pipe  and  observed  that  the  variation  of  the  Young's  modulus,  which  depends  on
temperature, plays a crucial role in our calculations. When a point load is applied with an imposed displacement, the Young's
modulus  does  not  significantly  affect  the  circumferential  deformation  because  the  resulting  deformation  is  primarily
controlled by the imposed displacement, regardless of the material's stiffness.

25
However, when internal pressure is added, the Young's modulus becomes important in changing the strain. Under internal
pressure,  the  material's  response  depends  on  its  elastic  properties,  particularly  the  Young's  modulus.  A  temperature-
dependent Young's modulus directly affects the material's stiffness and hence its ability to resist deformation under an applied
load. This means that different deformations will occur under the same pressure, depending on the temperature. Figure23.
illustrate how the changing of the young modules affect in this step.

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0.16

0.14

0.12

0.1

0.08

0.06

0.04

0.02

0

0

Article Referance
𝜶 Constant-1,5*10-4-NAllowSep-EV
𝜶 Constant-1,5*10-4-NAllowSep-EC

5

10

15

20

25

30

35

Inner fiber from the pin to the end [mm]

Figure23. The effect of the changes of young modulus in function de temperature in the Circumferential strain computed along the
inner fiber of the pipe wall starting below pin position after heating until 80°C and adding internal pressure.

Adding internal pressure highlights the importance of the Young's modulus, especially when it varies with temperature.
For materials subjected to simultaneous thermal and mechanical loads, it is essential to consider how their mechanical
properties, such as the Young's modulus, change with temperature. The results show that for accurate and realistic
simulations, it is crucial to model the Young's modulus as a function of temperature. Ignoring this variation can lead to
significant errors in predicting deformations and stresses under varied thermal and mechanical conditions.

Considering the variability of mechanical properties with temperature allows for better prediction of the behavior of
structures under complex loads, ensuring their reliability and safety.

Elasto-Plastic  Modeling: Using the tensile test data  for PE100,  we processed the information to perform our analysis in
Abaqus for each model, as illustrated in the corresponding figures. We considered three cases: perfect plasticity with a yield
strength of 27 MPa, linear plasticity modeled by selecting two points on the stress-strain curve, and a more realistic model
by intelligently choosing data points to approximate the material's actual behavior, Figuer24.

Perfect plasticity

Linear plasticity

Realistic model

140.0

120.0

100.0

80.0

60.0

40.0

20.0

0.0

0.0

140.0

120.0

100.0

80.0

60.0

40.0

20.0

0.5

1.0

1.5

2.0

0.0

140.0

120.0

100.0

80.0

60.0

40.0

20.0

0.0

0.0

0.5

1.0

1.5

2.0

0.0

0.5

1.0

1.5

2.0

Figure24. Different model elasticity that applied in simulation

The results of our simulations are presented below figure 25. There are three different curves, each corresponding to a specific
case. The curve corresponding to perfect plasticity with a threshold of 27 MPa shows a simple model to simulate but is less

26
realistic. The linear plasticity curve has a different shape from the other curves, with no peak. The real case curve, plotted
on this graph, shows that although the behavior is similar to that of perfect plasticity, the values are very different. Therefore,
it is preferable to choose the most realistic model to obtain results close to reality, even if the calculation time increases.

0.14

0.12

0.1

0.08

0.06

0.04

0.02

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0

0

5

Perfect plasticity-Yield  27

Linear  plasticity-Yield  27

Realistic model -Yield  27

10
25
Inner fiber from the pin to the end [mm]

15

20

30

35

Figure 25. The effect of the different hardening model in the Circumferential strain computed along the inner fiber of the pipe wall
starting below pin position after heating until 80°C

 We then compared our results with those from the reference article. Figurs26.  As observed, the curve with a 27 MPa yield
strength did not match the reference article's results. However, by using a higher yield strength of 35 MPa and a perfect
plasticity  model,  we obtained a curve that closely aligns  with the reference data  in both  the point load and point load +
heating scenarios. Despite this, it should be noted that this does not represent a real-world case.

0.14

0.12

0.1

0.08

0.06

0.04

0.02

)
-
(
n
i
a
r
t
S
l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0

0

article

Plasticité parfaite - LE 35

Perfect plasticity-Yield  27

5

10

15

20

25

30

Inner fiber from the pin to the end [mm]

0.16

0.14

0.12

0.1

0.08

0.06

0.04

0.02

)
-
(
n
a
r
t
S

i

l

a
i
t
n
e
r
e
f
m
u
c
r
i
C

0

0

15
5
Inner fiber from the pin to the end [mm]

20

10

25

30

Figurs26. Circumferential strain computed along the inner fiber of the pipe wall starting below pin position in PLT scenario ( left)
and PLT + Heating Scenario(Right)

 Subsequently, we introduced internal pressure and observed that with the adjusted yield strength, the results improved. A
parametric study was conducted by varying the yield strength and temperature (80°C). As depicted in the figure, increasing
the yield strength initially causes deformation to increase, and the peak of the curve shifts to the left, indicating a reduction
in the plastic zone. However, beyond a certain value, further increasing the yield strength at 80°C results in a decrease in
total deformation. Figure27.

27

Article

TD- LE35 to5 at 80°C

TD- LE35 to 15 at 80°C

TD- LE35 to 16 at 80°C

TD- LE35 to 20 at 80°C

)
-
(

n
i
a
r
t
S

l
a
i
t
n
e
r
e
f
m
u
c
r
i

C

0.2

0.18

0.16

0.14

0.12

0.1

0.08

0.06

0.04

0.02

0

0

5

10

15

20

25

30

35

Inner fiber from the pin to the end [mm]

Figure 27. Parametric Study on the Variation of Yield Strength at 80°C

4.2.  PHASE 2: CHEMICAL DEGRADATION ASSESSMENT: EFFECTS OF

CHLORINE DISINFECTION

In this section, we present the results of various tests conducted to assess the mechanical and physico-chemical properties of
HDPE, specifically PE100. The objective was to analyze these properties both before and after degradation due to chlorine
exposure. However, due to the extensive time required for aging on the platform provided by SUEZ, only the results before
aging will be presented at this time.

Additionally,  we  have  included  results  on  the  mechanical  and  physico-chemical  behavior  of  the  material  after  film
fabrication. The results are structured as follows:

Firstly, we present the mechanical behavior of HDPE using simple uniaxial traction tests in three states: the sample of the
pipe in its original form, samples fabricated from a film with a thickness of 600 microns at ambient temperature, and the
mechanical behavior of the film at 40°C.

Following  this,  we  will  detail  the  chemical  characterization  using  thermogravimetric  analysis  (TGA),  Fourier  transform
infrared spectrometry (FTIR), differential scanning calorimetry (DSC), and oxidative induction time (OIT).

Mechanical characterization: uniaxial tensile tests

Our study focused on determining the Young's modulus, elastic limit, and elongation at break of the materials. Tests were
conducted  under  controlled  conditions  of  23°C  and  50%  humidity,  using  different  test  speeds:  10  mm/min  for  thicker
specimens and 5 mm/min for films with a thickness less than 1 mm. For tests conducted at 40°C, the films were tested at the
same speed of 5 mm/min. Once the raw stress-strain curve is obtained, which shows the variation in load (F) applied to the
specimen as a function of the elongation (ΔL) it undergoes, it is necessary to normalize this curve to eliminate the influence
of the initial dimensions of the specimen and to focus solely on the material's properties.

To achieve this, we divide the load (F) by the initial cross-sectional area (S0) of the specimen to obtain the nominal stress
(σ). The vertical axis is then scaled in stress units (MPa). Similarly, we divide the elongation (ΔL), which is the difference
between the instantaneous length (L) and the initial length (L0), by the initial length (L0). This gives us the strain, typically
expressed as a percentage (%), and the horizontal axis becomes the strain axis (ε). This normalized curve is referred to as the
nominal tensile stress-strain curve or simply the tensile curve.

However, for inputting material properties into Abaqus, we require true stress and true strain curves. Therefore, further data
processing is necessary. To obtain the true stress (𝜎𝑡𝑟𝑢𝑒) and true strain (𝜀𝑡𝑟𝑢𝑒) required for more accurate analysis, we use
the following transformations:

𝜎𝑡𝑟𝑢𝑒 =  𝜎 (1 + 𝜀)

𝜀𝑡𝑟𝑢𝑒 =  𝐿𝑛 (1 + 𝜀)

Eq.9

Eq.10

28
These calculations ensure that the stress-strain data accounts for the continuous deformation and actual material behavior
during testing, providing a more precise representation of the material's response under load.

In Figure 28, I will present the results obtained for PE100 in both pipe and film forms with varying thicknesses. It is evident
that  samples  fabricated  from  the  film  using  the  compression  press  method  at  180°C  and  160°C  exhibit  nearly  identical
mechanical  properties.  This  observation  allows  us  to  conclude  that  reheating  the  material  to  a  molten  state  and  then
fabricating the film does not significantly alter its mechanical properties.

Figure 28. Results of the Tensile Test for PE100: Nominal Stress-Strain Curve (left) and True Stress-Strain Curve (right)

In the figure29. you can observe the properties of PE100 at 40°C. However, due to equipment limitations within the oven,
measurements for elongation at yield and yield strength were not possible. Nevertheless,  we have captured the material's
behavior within the plastic range, which is essential for inputting data into Abaqus.

In the table.5. all the mechanical properties information is shown. Some information for the test at 40 °C is missing due to
machine limitations (limited imposed displacement).

Figure 29. Results of the Tensile Test for PE100 at 40°C: Nominal Stress-Strain Curve (left) and True Stress-Strain Curve (right)

Exposure time: 0 Day
Semple

PE100-Pipe -001
PE100-Pipe -002

Young's Modulus
(MPa)
751,3
668,0

Elongation at Break
(%)
861,4
841,3

Yield Strength
(MPa)
21,45
20,40

Tensile Strength
(MPa)
26,03
25,58

PE100-Pipe -003
PE100-Film-180°C-001 (at 23°C)
PE100-Film-180°C -002 (at 23°C)
PE100-Film-180°C -003 (at 40°C)
PE100-Film-180°C -004 (at 40°C)
PE100-Film-180°C -005 (at 40°C)
PE100-Film-180°C -006 (at 40°C)
PE100-Film-180°C -007 (at 40°C)

ATG Results

757,2
782,4
772,6
511,5
490,8
460,9
331,3
489,4

890,9
1004,8
1493,0

21,27
20,00
20,60
12,1
13,3
15,3
15,3
14,2

Table5. the mechanical properties of PE100

29

27,28
33,11
30,61

Thermogravimetric  analysis  (TGA)  is  a  thermal  technique  well-suited  for  assessing  the  stability  of  polymeric  materials.
Under high temperatures, polymers decompose, resulting in the formation of various low molecular weight products. TGA
provides valuable insights into the behavior of polymers during oxidation or thermal degradation, offering predictive data
for  their  performance  under  real  atmospheric  conditions.  In  this  experiment,  we  anticipate  that  TGA  results  will  reveal
structural changes due to degradation from chlorine disinfectant exposure.

The thermal stability of PE100 composites was initially investigated using TGA under a nitrogen atmosphere with a heating
rate of 20°C/min. The analysis revealed that the decomposition of HDPE-PE100 composites occurs in two distinct steps
(Table 6). The first decomposition step occurs approximately between 485°C and 490°C, corresponding to the degradation
of Polymer part. The second degradation step, around 615°C to 625°C, is attributed to the  degradation of HDPE carbon
black. The residual mass after polymer degradation was found to be approximately 1.8% of the initial mass, indicating that
the  carbon  black  content  in  PE100  is  1.8%  of  the  initial  mass,  while  the  remaining  98.2%  consists  of  polymer  and
antioxidants. Table 6 summarizes the temperatures of polymer and carbon black degradation, along with the corresponding
mass percentages before and after film fabrication of PE100, providing crucial data for the ATG section of my report.

Exposure time: 0 Day

specimens

Temp _Deg of polymer (°C)
Temp _Deg of carbon black (°C)

PE100-Film-180°C-
001
491,9
616,6

PE100-Film-160°C-
001
489,7
624,2

PE100-Film-160°C-
002
490,3
622,3

% Mass of Polymer

% Mass of carbon black

98, 2%

1,8%

98,2%

1,8%

98,1%

1,9%

Table 6. the result of the temperatures of polymer and carbon black degradation

PE100-Pipe-001

485,8
620,6

98,2%

1,8%

FTIR spectroscopy results

In this section, I will present the results of the experimental characterization, specifically focusing on the FTIR spectroscopy
results. As previously explained, we characterized the materials at different stages. Initially, we characterized the materials
in their non-aging and virgin stages (PE100). Once we manufactured the films under high temperature (180°C) and pressure,
we performed another round of characterization to ensure that these conditions did not alter the properties of the materials.
Therefore, all the data presented here are before and after the fabrication of...

Figures30. depict the transmittance spectra of infrared waves for the samples before and after film production. These spectra
will be analyzed to interpret any changes in the chemical structure of the materials throughout the fabrication process.

The observed peaks at 2914, 2847, 1743, 1461, and 716 cm⁻¹ can be interpreted as follows:

•  Peaks at 2914 cm⁻¹ and 2847 cm⁻¹ indicate C-H stretching vibrations, confirming the presence of methylene groups

characteristic of polyethylene.

•  A  peak  at  1743  cm⁻¹  suggests  the  presence  of  carbonyl  groups,  indicative  of  some  oxidation,  possibly  due  to

environmental exposure.

•  The peak at 1461 cm⁻¹ corresponds to CH2 bending vibrations, further confirming the polyethylene structure.
•  The peak at 716 cm⁻¹ is associated with CH2 rocking vibrations, indicating the crystalline nature of the polymer.

30

2914

2847

PE100-Pipe

PE100-Film160°C

PE100-Film180°C

716

1461

1743

500

1000

1500

2000

2500

3000

3500

4000

4500

cm-1

Figures30. FTIR spectroscopy results for PE100 before and after film fabrication

0.45

0.4

0.35

0.3

0.25

0.2

0.15

0.1

0.05

0

0

DSC Results

DSC measurements were conducted to monitor changes in the crystallinity of HDPE pipes before and after film fabrication
at  temperatures  of  180°C  and  160°C.  This  was  important  due  to  the  known  effect  that  structural  changes,  induced  by
fabrication or degradation, can have on material properties. From the results presented in Table 7, it is evident that the melting
temperatures and crystallinity rates were recorded.

Across all samples, there was no significant variation in melting temperatures. However, a slight difference in crystallinity
rates was observed specifically in the samples fabricated at 160°C. It is hypothesized that the rapid cooling method used may
have influenced these crystallinity rates.

Exposure time: 0 Day
PE100
T_Melting

crystallinity rates
PE100- Film-180°C -Slow cooling
T_Melting
crystallinity rates
PE100- Film-160°C-Rapid Cooling
T_Melting
crystallinity rates
PE100- Film-160°C - Slow cooling
T_Melting
crystallinity rates

Specimens 1
129,8

Specimens 2
129,9

63,2%

64,6%

Specimens 3
129,3

64,4%

129,6
61,4%

129,5
61,2%

127,9
57,2%

128.7
61.4

129,8
60,7%

128,3
57,3%

129.7
63.3

Moyenne
129,7

64,1%

129,7
61,1%

128,1
57,3%

129.2
62.35

Conversely, for the PE films fabricated at 160°C with a slow cooling process, the crystallinity rate increased and closely
matched that of the original PE100 pipe before film fabrication.

Table 7. Result of DSC test for PE100 in different condition

OIT Results

This section presents the results of the Oxidative Induction Time (OIT) tests conducted using two different methods. In the
first method, the material was heated to 220°C, and in the second method, it was heated to 200°C. All results are summarized

31
in Table 8. As shown, the OIT values increase as the temperature decreases. This observation aligns with the Arrhenius law,
which states that the reaction rate decreases with a reduction in temperature. Consequently, the rate of oxidation slows down
at lower temperatures, resulting in higher OIT values.

these results indicate that the oxidative stability of the material improves at lower temperatures. The increase in OIT with
decreasing temperature suggests that the material is less prone to oxidation and degradation under lower thermal conditions.
This behavior is crucial for applications where long-term thermal stability is required, such as in piping systems exposed to
varying environmental conditions.

The variation in oxidation rates with temperature also provides valuable insights for predicting the material's performance
and lifespan. By understanding the temperature dependence of oxidation, we can better assess the durability and reliability
of HDPE pipes and films in real-world applications.

Exposure time: 0 Day
Samples

PE100-Film160°C-001
PE100-Film180°C-001
PE100-Pipe-001
PE100-Pipe-002
PE100-Pipe-003
PE100-Pipe-004
PE100-Pipe-005

mass
(mg)
5.07
5.39
7.94
6.24
8.88
8.16
7.35

Test temperature
(°C)
220
220
220
220
200
200
200

OIT
(min)
21.3
19.5
22.4
23.0
122.7
121.0
116.5

The oxidation rate
(W/g.min)
1.
0.7
0.7
0.9
0.1
0.1
0.1

Table8. The results of OIT for the PE100 in different condition

5. Concluding Remarks
The research undertaken during this internship involved a comprehensive analysis of the performance and lifetime of PE100
pipes, particularly under the combined effects of mechanical punching and chemical degradation. Before this, it was crucial
to validate our modeling with a reference case to ensure that the chosen parameters were optimal among the many options
available in Abaqus. Therefore, the first step was to validate our model by simulating a case from an existing article and then
comparing the results with those in the article.

In this purpose, Abaqus Finite Element software was utilized to simulate the mechanical stresses and strains in PE100 pipes
under  various  conditions.  A  quarter  of  a 3D  model  of  the  pipe, the  pin  and  the  counter  support  were  simulated  in  three
sequence  steps  (Punching  + Heating  + Pressure intern).  In this  part of the  project, we  understood that  when  an imposed
displacement is applied, the value of Young's modulus is not as crucial as when an imposed force is applied. Additionally,
in the case of point load, it is not feasible to reduce our model to a 2D case due to the inability to capture the complex three-
dimensional stress distribution around the punch. Significant out-of-plane stresses and deformations arise during punching,
which a 2D model cannot accurately represent. In contrast, the 3D model effectively captures these complexities, resulting
in better agreement with the reference results.

It was also observed that thermal deformation can be as significant  as mechanical  deformation when  the pipe operates at
high temperatures. In the case of punching load and heating, without considering the coefficient of thermal expansion, the
results  remain  unchanged.  However,  when  including  the  change  in  Young's  modulus  as  a  function  of  temperature,  the
deformation under the pin shows significant variation. The initial value at zero point is 0.81, and after adding the coefficient
of thermal expansion, the value changes to 0.1, indicating that the deformation altered by more than 20%.

In the case of elasto-plastic behavior, it was also observed that when adding the plasticity model, the change in deformation
is not linear as in the elastic case. The variation in deformation is significantly different depending on the yield limit value,
not only at ambient temperature but also at 80°C.

The study  compared three plasticity  models:  perfect plasticity,  linear plasticity,  and a model  approximating  real material
behavior.  It  was  evident  that  the  perfect  plasticity  model,  despite  its  simplicity,  fails  to  capture  the  nuances  of  material
deformation  under  stress.  The  linear  plasticity  model,  although  closer  to  reality,  still  diverged  significantly  from  actual

32
behavior.  The  most  realistic  model,  which  involves  strategically  selecting  data  points,  provided  the  most  accurate
representation, albeit with higher computational costs.

Additionally,  in the case of PL + heating  + internal  pressure, it was  observed that when  adding the plasticity  model, the
change in deformation is not linear as in the elastic case. The deformation varies significantly depending on the yield limit
value, not only at ambient temperature but also at 80°C.

In the second phase of the project, we would investigate the chemical degradation of polymer pipes exposed to chlorinated
disinfectants commonly used in drinking water treatment. Based on articles, we predicte that changes in the number of chains
and, consequently, molecular weight caused by chain scissions can significantly alter the physical and mechanical properties
of PE pipes. However, since we do not yet have the aging materials for the tests, we cannot present any results in this report.
The aging tests were initiated on May 26, and we do not yet have the specimens after one month of aging to observe this
effect.

Thus, we only have results for the virgin material. In the investigation of the mechanical properties of PE100, it was observed
that the properties of the material for the film  fabricated at 180°C or at 160°C do not change significantly.  However, the
cooling rate in the process of fabrication of the film can impact the crystallinity of the material. Rapid cooling results in a
slightly lower crystallinity rate. therefore, consequently its mechanical properties will change.

With doing the tensile test at tow temperature we observed that, the properties of the material at 40°C during the tensile test
change drastically. For example, the Young's modulus at 23°C is on average 746 MPa, while at 40°C it drops to 456 MPa,
representing a decrease of up to 40%. The yield strength also decreases from 21 MPa to 14 MPa, a change of 32%.

Finally, the section concludes with a consideration of lifetime prediction, specifically focusing on the combined influence of
punching and chemical degradation in polymer pipes. The literature review highlights the gaps in assessing the lifetime of
polymer pipes, specifically HDPE, concerning the combination of two phenomena: punching caused by external load and
chemical decay due to chlorine disinfectants in water systems.

As  the third phase of the project depends on the  data from  phase two  and the results  after degradation, we  are currently
unable to proceed with this phase

Moreover, while our comparisons thus far have been primarily visual, it is imperative to enhance precision and reliability
through  analytical  comparisons  in  future  studies.  Refining  computational  models  to  integrate  intricate  degradation
mechanisms and their interactions is essential.

Our  study  proposes  the  use  of  computational  models  and  chemical  assessments  for  a  detailed  analysis.  The  goal  is  to
understand the mechanical  stresses from different PL parameters and estimate  pipe lifetime, considering both mechanical
and chemical degradation. This approach is crucial for developing predictive models and maintenance strategies, ultimately
extending the lifetime of polymer pipes in water distribution systems.

References

33

A. Boujlal, Analysis Of Slow Crack Initiation Of Old Polyethylene Resins By Means Of An Elasto-Visco-Plastic
Rheological Model: Experimental And Numerical Approach, Plastic Pipes XVI, Barcelona, Spain, 2012

A. Frank, “A numerical methodology for lifetime estimation of HDPE pressure pipes”, Engineering Fracture

Mechanics, 78, 2011, p. 3049–3058.

A. Frank, G. Pinter, R.W. Lang, “Prediction of the remaining lifetime of polyethylene pipes after up to 30 years in

use”, Polymer Testing 28, 2009, p. 737–745.

A. Sukhadia, “Assessing the Slow Crack Growth Resistance of PE Resins and Pipe Service Lifetimes Predictions”,

15th Plastic Pipes conference, PPXV, Vancouver, Canada, 2010.

A. Tripath, S. Mantell, J. Le, “Chemo-mechanical modeling of static fatigue of high density polyethylene in bleach

solution”, International Journal of Solids and Structures, 217–218, 202, p. 90–105.

ASTM D2837, “Standard Test Method for Obtaining Hydrostatic Design Basis for Thermoplastic Pipe Materials or
Pressure Design Basis for Thermoplastic Pipe Products”, ASTM D2837, ASTM International, West
Conshohocken, 2013.

B.-H. Choi, A. Chudnovsky, R. Paradkar, W. Michie, Z. Zhou, P.-M. Cham, "Surface degradation of HDPE-100 pipe:

Effects of some aggressive environments (solvents) ", Polymer Degradation and Stability, 94 (5), 2009, p.
859–867 May 2009.

B. Fayolle, X. Colin, L. Audouin, J. Verdu, “Mechanism of degradation induced embrittlement in polyethylene”,

Science Direct Polymer Degradation and Stability, 92, 2007, p. 231-238.

DIN EN ISO 9080, “Plastics piping and ducting systems –Determination of the long-term hydrostatic strength of

thermoplastics materials in pipe form by extrapolation”, 2003.

Dossier technique conduit 9010 RC en PE 100-RC.

E. Van Der Stok, F. Scholten, “Two complementary tests to determine the quality of pe 100-RC”, 17th Plastic Pipes

Conference PPXVII, Chicago, Illinois, USA, 2014.

E.M. Hoang, D. Lowe, “Lifetime prediction of a blue PE100 water pipe”, Polymer Degradation and Stability, 93(8),

2008, p. 1496.

F.L. Scholten, D. Gueugnaut, F. Berthier, “A More Reliable Detergent for Cone and Full Notch Creep Testing of PE

Pipe Materials”, 11th Plastic Pipes conferences, PPXI, Munich,Germany, 2001.

H. Bromstrup, PE 100 pipe systems, 2nd edition, 2004, p. 15- 35.

ISO 13480, “Polyethylene pipes - Resistance to slow crack growth - Cone test method”, 1997.

ISO 9080, Plastics Piping and Ducting Systems, “Determination of the Long-term Hydrostatic Strength of

Thermoplastics Materials in Pipe Form by Extrapolation”, ISO 9080, CEN, Brussels, 2003.

J. Hassinen, M. Lundbäck, M. Ifwarson, U.W. Gedde, “Deterioration of polyethylene pipes exposed to chlorinated

water”, Polymer Degradation Stability, 84(2), 2004, p.261-267.

J. Hessel, “Minimum service-life of buried polyethylene pipes without sand embedding”, 3R international Special

Plastics Pipes, 40(6), 2001, p. 178-184.

J. Lenz, “Long-Term Point Loads Testing with Plastic Pipes”, 11th plastic pipe conference, PPXI, Munich, Germany,

2001.

L. Andena, M. Rink, R. Frassine, and R. Corrieri, “A fracture mechanics approach for the prediction of the failure time

of polybutene pipes”, Engineering Fracture Mechanics, 76(18), 2009, p. 2666-2677.

M. Farshad, “Two new criteria for the service life prediction of plastics pipes”, Polymer Testing, 23(8), 2004, p. 967-

972.

M. Haager, “Fracture mechanics methods for the accelerated characterization of the slow crack growth behavior of
polyethylene pipe materials”, Doctoral Dissertation, Institute of Materials Science and Testing of Plastics,
University of Leoben, Austria, 2006.

34

N. Brown, X. Lu, “Controlling the Quality of PE Gas Piping Systems by Controlling the Quality of the Resin”,
Proceedings of the 13th Plastic Fuel Gas Pipe Symposium, San Antonio, Texas, USA, 1993.

N. Brown, X. Lu, “PENT quality control test for PE gas pipes and resins”, Proceedings of the 12th Plastic Fuel Gas

Pipe Symposium, Boston, Massachusetts, USA, 1991.

P. Hutař, M. Ševčík, L. Náhlík, G. Pinter, G. Pluvinage, M. Elwany, “Safety, Reliability and Risk Associated with

Water, Oil and Gas Pipeline”, 2007.

P. Hutar, M. Sevcik, L. Nahlik, G. Pinter, A. Frank, I. Mitev, “A numerical methodology for lifetime estimation of

HDPE pressure pipes”, Engineering Fracture Mechanics, 78(17), 2011, p. 3049-3058P.

P. Hutar, M. Ševcík, A. Frank, L. Náhlík, J. Kucera, and G. Pinter, “The Effect of Specimen Size on the Determination
of Residual Stress in Polymer Pipe Wall”, Engineering Fracture Mechanics, 108, 2013, p. 98.PAS1075:2009-
04, “Pipes made from Polyethylene for alternative installation techniques – Dimensions, technical
requirements and testing”, 2009-04.

Plastics Pipe Institute (PPI), Handbook of Polyethylene Pipe, 2nd edition,2008, p. 43-103.

S. Nestelberger, J. Cheng, “Finite element method (fem) used to simulate the stress/strain of the point load test (PLT)
in a broader study to support the future iso test standard”, Proceedings of the 20th Plastic Pipes Conference
PPXX, Amsterdam, Netherland, 2021.

U. Niebergall, O. Mertlová, E. Nezbedová, “Full Notch Creep Test - ISO Round Robin Test”, 13th Plastic Pipes

conference, PPXIII, Washington DC, USA, 2006.

V. Rouyer, M. Cornette, “Resistance of Crosslinked PE Pipes to Rock Impingement”, Proceedings of the 11th Plastic

Pipes Conference PPXI, Munich Germany, 2001.

W. Yu, B. Azhdar, D. Andersson, T. Reitberger, J. Hassinen, T. Hjertberg, U.W. Gedde, “Deterioration of

polyethylene pipes exposed to water containing chlorine dioxide”, Polymer Degradation and Stability, 96,
2011, p. 790-797.

X. Colin, L. Audouin, J. Verdu, M. Rozental-Evesque, B. Rabaud, F. Martin, F. Bourgine, “Aging of Polyethylene

Pipes Transporting Drinking Water Disinfected by Chlorine Dioxide. I. Chemical Aspects”, Polymer
Engineering And Science, 49(7), 2009, p. 1429-14376.

X. He, X. Zha, X. Zhu, X. Qi, B. Liu, “Effect of short chain branches distribution on fracture behavior of polyethylene

pipe resins”, Polymer Testing, 68, 2018, p. 219-228.

X. Zheng, X. Zhang, L. Ma, W. Wang, J.g Yu, “Mechanical characterization of notched high density polyethylene
(HDPE) pipe: Testing and prediction”, International Journal of Pressure Vessels and Piping,173, 2019, p.
11-19.

Y. Huang, Q. Zhang, X. Lu, Y. Gong, H. Zhou, J. Feng, “Comparative Investigation on Step-cycle Tensile Behaviors

of Two Bimodal Pipe-grade Polyethylene with Different Slow Crack Growth Resistance”, Chinese Journal of
Polymer Science, 38, 2020, p. 611–619.

Z. J. Zho, D. Chang, “Stress intensity effect on slow crack growth in scratched polyethylene pipe”, 15th Plastic Pipes

Conference PPXV, Vancouver, Canada, 2010.

Z. Xiong, R. Wei, “A New Point Loading Test Method For HDPE Pipe Materials”, 16th Plastic Pipes Conference,

PPXVI, Barcelona, Spain, 2012.


