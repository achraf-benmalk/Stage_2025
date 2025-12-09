Modélisation de la Dégradation et
Prédiction de la Durée de Vie de
Canalisations en Polyéthylène Exposées
aux Désinfectants Chlorés

Stage de Master 2 Recherche
Spécialité : Mécanique des Matériaux pour l’Ingénierie et l’Intégrité des

Structures (MAGIS)

Réalisé par :

Benmalk Achraf

Encadré par :

Xavier COLIN
Juan-Pablo MARQUEZ-
COSTA

(Arts et Métiers – Laboratoire PIMM)
(Arts et Métiers – Laboratoire PIMM)

Année universitaire 2024-2025

Résumé

Ce stage avait pour objectif d’évaluer la pertinence d’un modèle cinétique de dé-
gradation existant, développé pour le système PE80/DOC, dans le contexte industriel
d’un grade PE100 exposé à l’hypochlorite de sodium (HOCl). Pour ce faire, le modèle
a été implémenté numériquement en Python, validé par rapport aux données de la
littérature, puis confronté aux résultats de caractérisations expérimentales, menées au
laboratoire sur des films vieillis fournis par SUEZ. L’implémentation numérique s’est
avérée fonctionnelle et capable de reproduire avec une bonne fidélité les résultats de
référence, validant ainsi notre outil de simulation. La confrontation avec les données
de vieillissement en eau pure a permis de confirmer que le modèle, avec ses paramètres
de transport d’origine, capture correctement la cinétique de perte physique d’antioxy-
dant. Pour le vieillissement en eau de Javel, une approche exploratoire utilisant une
concentration "DOC effective" a montré que le modèle pouvait simuler qualitativement
l’accélération de la dégradation. Cependant, les écarts quantitatifs et le comportement
non-monotone de l’indice carbonyle expérimental soulignent les limites de cette ana-
logie et la nécessité de développer un schéma cinétique spécifique à l’HOCl. Enfin,
l’analyse des propriétés mécaniques a mis en évidence une corrélation expérimentale
directe et très nette entre la dégradation chimique (chute de l’OIT) et la fragilisation
du matériau (chute de l’allongement à la rupture). Ce résultat clé valide qualitative-
ment l’hypothèse fondamentale du couplage mécano-chimique sur laquelle repose le
modèle et ouvre des perspectives prometteuses pour la prédiction de durée de vie.

Abstract

This internship aimed to assess the relevance of an existing kinetic degradation
model, originally developed for the PE80/DOC system, in the industrial context of
a PE100 grade exposed to sodium hypochlorite (HOCl). To achieve this, the model
was numerically implemented in Python, validated against literature data, and then
confronted with experimental characterization results, carried out in the laboratory on
aged films provided by SUEZ. The numerical implementation proved to be functional
and capable of faithfully reproducing the reference results, thereby validating our si-
mulation tool. The confrontation with aging data in pure water confirmed that the mo-
del, with its original transport parameters, correctly captures the kinetics of physical
antioxidant loss. For bleach aging, an exploratory approach using an "effective DOC"
concentration showed that the model could qualitatively simulate accelerated degra-
dation. However, quantitative discrepancies and the non-monotonic behavior of the
experimental carbonyl index highlight the limitations of this analogy and the need to
develop a specific kinetic scheme for HOCl. Finally, the analysis of mechanical proper-
ties revealed a direct and clear experimental correlation between chemical degradation
(OIT drop) and material embrittlement (drop in elongation at break). This key result
qualitatively validates the fundamental hypothesis of the chemo-mechanical coupling
on which the model is based and opens promising prospects for lifetime prediction.

Table des matières

Résumé

Abstract

1

Introduction : Motivation de l’étude

2 Étude Bibliographique du Problème

2.1 Résultats Acquis : Mécanismes et Modèles de Dégradation . . . . . . . .
2.1.1 Rupture Différée des Polyéthylènes sous Contrainte . . . . . . . .
. . . . . . . .
2.1.2 Le Modèle Cinétique Multi-Niveaux de Colin et al.
. . . . . . . . .
2.1.3 Modèles Chemo-Mécaniques pour l’Eau de Javel
2.2 Points Durs et Piste de Recherche . . . . . . . . . . . . . . . . . . . . . . .
. . . . . . . . . . . . . . . . . . . . . . .
2.2.1 Les Verrous Scientifiques
2.2.2 Piste de Recherche et Objectifs du Stage . . . . . . . . . . . . . . .

3 Démarche Scientifique et Outils Numériques

3.2

3.1 Formalisation Mathématique du Modèle . . . . . . . . . . . . . . . . . . .
Système d’Équations de Réaction-Diffusion . . . . . . . . . . . . .

3.1.1
3.1.2 Description du Problème Physique et Conditions aux Limites (CL)
Implémentation et Résolution Numérique . . . . . . . . . . . . . . . . . .
3.2.1 Discrétisation Spatiale : la Méthode des Lignes . . . . . . . . . . .
3.2.2
Intégration Temporelle et Gestion de la Raideur . . . . . . . . . .
3.2.3 Calcul des Grandeurs d’Intérêt . . . . . . . . . . . . . . . . . . . .
3.3 Matériaux et Méthodes Expérimentales . . . . . . . . . . . . . . . . . . .
3.3.1 Matériau et Préparation des Échantillons . . . . . . . . . . . . . .
3.3.2 Conditions de Vieillissement
. . . . . . . . . . . . . . . . . . . . .
3.3.3 Techniques de Caractérisation . . . . . . . . . . . . . . . . . . . .
Stratégie de Confrontation Modèle/Expérience . . . . . . . . . . . . . . .

3.4

4 Résultats et Analyses

4.1

Stratégie de Modélisation et Paramètres de Simulation . . . . . . . . . .
. . . . . . . . .
4.1.1 Paramètres du Système de Référence (Validation)

3

1

2

1

2

2
2
3
6
6
6
6

8

8
8
9
9
9
10
12
12
12
13
13
14

16

16
16

17
18
18

18
19
19

20

22

24

26

lin et al. 2009) .

.
4.4 Confrontation avec les Données SUEZ (Films PE100)

4.2 Analyse des Profils de Concentration Simulés (Résultats Préliminaires) .
4.3 Validation du Modèle Numérique par Confrontation à la Littérature . .
4.3.1 Validation des Profils d’OIT (Réf. Figure 3, Colin et al. 2009)
. . .
4.3.2 Validation des Profils de Carbonyles et d’OIT (Réf. Figure 4, Co-
. . . . . . . . . . . . . . . . . . . . . . . . . . . .
. . . . . . . . . . .
4.4.1 Vieillissement en Eau Pure (H2O) : Cinétique de Perte Physique .
4.4.2 Vieillissement en Eau de Javel (HOCl) : Dégradation Chimique
. . . . . . . . . . . . . . . . . . . . . . . . . . . .
Accélérée
Impact sur les Propriétés Mécaniques : Vers une Prédiction de
. . . . . . . . . . . . . . . . . . . . . . . . . . . .
Durée de Vie

4.4.3

.

.

.

.

.

.

5 Conclusions et Perspectives

Bibliographie

4

Table des figures

après exposition au DOC. Adapté de Yu et al. [8].

2.1 Régimes de rupture du PE . . . . . . . . . . . . . . . . . . . . . . . . . . .
2.2 Micrographie d’une fissure fragile initiée en paroi interne d’un tuyau PE
. . . . . . . . . . . . .
2.3 Accumulation de produits chlorés ([PCl]) en fonction du temps d’expo-
sition dans une solution de dioxyde de chlore, illustrant la réaction de
. . . . . . . . . . . .
greffage du désinfectant. Adapté de Colin et al. [2].
2.4 Formation de produits carbonylés ([CO]), signature de l’oxydation de
la matrice, dans un film de PE exposé au dioxyde de chlore. Adapté de
. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Colin et al. [2].
Illustration de la transition ductile-fragile : chute de l’allongement à la
rupture (εr) lorsque la masse molaire (Mw) passe en dessous d’un seuil
. . . . . . . . . . . . . . . . . . . . .
critique. Adapté de Fayolle et al. [5].
2.6 Simulation de courbes contrainte-durée de vie pour du PEHD, montrant
l’accélération de la rupture en présence d’eau de Javel. Adapté de Tripa-
. . . . . . . . . . . . . . . . . . . . . . . . . . . . .
thi et al. [7].

2.5

. .

. .

.

.

.

.

.

.

.

3.1

Illustration du principe d’une méthode de Runge-Kutta d’ordre 2 (RK2),
un schéma explicite performant pour les problèmes non-raides. Adapté
du support de cours MASC-FEM [6]. . . . . . . . . . . . . . . . . . . . . .
3.2 Approche globale du projet, montrant l’interconnexion entre les don-
nées d’entrée (A), les blocs de simulation (B), les prédictions (C) et les
. . . . . . . . . . . . . . . .
données de référence pour la validation (D).

scénarios de vieillissement à 40◦C.

4.1 Profils de concentration simulés pour l’Oxygène, l’Antioxydant (AH) et
. . . .
le Dioxyde de Chlore (DOC) à différents temps de vieillissement.
4.2 Confrontation des profils d’OIT simulés (a) et expérimentaux (b) pour 5
. . . . . . . . . . . . . . . . . . . . . .
4.3 Confrontation des profils d’OIT (pointillés) et de Carbonyles (symboles
. . . . . .
4.4 Évolution de l’OIT vs. Temps pour films PE100 en eau pure à 40 ◦C.
. . . . . . .

pleins) simulés (a) et expérimentaux (b) après 99 jours à 40◦C.

Confrontation : Données SUEZ (bleu) et Simulation (rouge).

5

2

3

3

4

5

6

11

14

17

18

19

19

4.5 Évolution de l’OIT (axe gauche) et de l’Indice Carbonyle (axe droit) en
fonction du temps pour des films de PE100 sous HOCl à 40 ◦C. Données
SUEZ : OIT (vert), Indice Carbonyle (orange). Simulation : OIT (violet),
[CO] en surface z=0 (marron). . . . . . . . . . . . . . . . . . . . . . . . . .
4.6 Évolution de la bande d’absorption carbonyle ( 1720 cm−1) dans les
spectres FTIR des films de PE100 en fonction du temps de vieillissement
. . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
sous HOCl. . .
4.7 Corrélation entre la dégradation chimique (OIT, axes de droite) et la dé-
gradation mécanique (Allongement à la rupture, axes de gauche) pour
des films de PE100 vieillis à 40◦C en eau pure et en eau de Javel (HOCl).

. .

. .

20

21

22

6

Chapitre 1

Introduction : Motivation de l’étude

Le polyéthylène (PE) est le matériau prédominant pour les réseaux de distribution
d’eau potable, apprécié pour sa durabilité, sa flexibilité et sa résistance à la corrosion.
Les grades modernes comme le PE100 sont conçus pour une durée de vie de plus de 50
ans, un enjeu économique et sanitaire majeur pour les opérateurs de réseaux comme
SUEZ.

Cependant, la garantie de cette longévité est remise en question par la nécessité
de désinfecter l’eau. L’utilisation d’agents chlorés, comme le dioxyde de chlore (DOC)
ou l’hypochlorite de sodium (HOCl), a été corrélée à des défaillances prématurées de
canalisations. Ces ruptures, souvent de nature fragile, suggèrent qu’un mécanisme de
dégradation chimique, initialement non pris en compte dans les modèles de durée de
vie purement mécaniques, joue en réalité un rôle critique.

Le verrou scientifique réside donc dans la capacité à développer et valider des mo-
dèles de prédiction de durée de vie qui intègrent ces phénomènes de vieillissement
chimique. Le laboratoire PIMM a développé un modèle mécano-chimique de référence
pour le système PE80/DOC [1]. Ce modèle constitue une base théorique solide, mais
sa validité pour d’autres systèmes reste une question ouverte.

L’objectif de ce stage est de relever ce défi en confrontant ce modèle à un double en-
jeu de transposition : évaluer sa pertinence et d’identifier ses limites pour un système
différent, un grade plus récent, le PE100, exposé à un désinfectant très répandu, l’hy-
pochlorite de sodium (HOCl). La démarche consiste à implémenter numériquement le
modèle en Python et à le confronter à un jeu de données de vieillissement expérimen-
tales inédites fournies par SUEZ.

1

Chapitre 2

Étude Bibliographique du Problème

2.1 Résultats Acquis : Mécanismes et Modèles de Dégra-

dation

2.1.1 Rupture Différée des Polyéthylènes sous Contrainte

La durée de vie des canalisations en PE est traditionnellement évaluée par des es-
sais de rupture sous pression hydrostatique constante. Les résultats sont présentés sous
forme de courbes de régression (Figure 2.1), qui révèlent typiquement trois régimes de
rupture [3].

FIGURE 2.1 – Schéma des courbes contrainte-durée de vie pour le PE, illustrant les ré-
gimes de rupture : (I) Ductile, (II) Fragile différée (mécanique), et (III) Fragile accélérée
(chimique).

Le Régime II, ou rupture fragile "physique", est gouverné par la propagation lente
de fissure (Slow Crack Growth, SCG). Ce mécanisme complexe, initié à partir de dé-
fauts, implique la cavitation de la phase amorphe et l’étirage des macromolécules jus-
qu’à la rupture des chaînes d’attache reliant les cristallites [5]. En l’absence d’agression
chimique, ce régime définit une durée de vie supérieure à 50 ans.

En présence de désinfectants, un Régime III apparaît, caractérisé par une chute
verticale de la durée de vie. Cette rupture prématurée est la conséquence directe de la
dégradation chimique du polymère (Figure 2.2).

2

FIGURE 2.2 – Micrographie d’une fissure fragile initiée en paroi interne d’un tuyau PE
après exposition au DOC. Adapté de Yu et al. [8].

2.1.2 Le Modèle Cinétique Multi-Niveaux de Colin et al.

Pour prédire l’apparition du régime III, le modèle de Colin et al. [1] couple la chimie
et la mécanique via une approche multi-niveaux, où chaque niveau décrit un phéno-
mène à une échelle différente.

Niveau 1 : Schéma Réactionnel et Produits d’Oxydation

Le cœur du modèle est un schéma cinétique d’autoxydation radicalaire en chaîne,
qui décrit la formation et la consommation des espèces chimiques clés. L’oxydation est
initiée par le désinfectant (DOC) et se propage via des radicaux peroxyles (P O•
2). Ce
processus génère divers produits d’oxydation, notamment des produits chlorés ([PCl])
par greffage direct du désinfectant sur les chaînes polymères (réaction k4d), comme
l’illustre expérimentalement la Figure 2.3.

FIGURE 2.3 – Accumulation de produits chlorés ([PCl]) en fonction du temps d’ex-
position dans une solution de dioxyde de chlore, illustrant la réaction de greffage du
désinfectant. Adapté de Colin et al. [2].

Le schéma réactionnel complet, incluant l’initiation, la propagation, la ramification

3

et la terminaison, est présenté ci-dessous :

Initiation : Y • + PH

k1d−−→ P • + Produits inactifs

Ramification : POOH

2POOH

2 + γ1sS + γ1coCO + H2O

k1u−−→ 2P • + γ1sS + γ1coCO
k1b−→ P • + P O•
k2−→ P O•
2

Propagation : P • + O2

P O•

2 + PH

k3−→ POOH + P •

Terminaison : P • + P • k4−→ Produits inactifs + γ4X

P • + P O•
2

P O•

2 + P O•
2

k5−→ Produits inactifs + γ5X
k6−→ Produits inactifs + O2
k7−→ POOH + A•

Piégeage AO : P O•

2 + AH

Y • + AH

k8d−−→ Y H + A•

Autres : P • + Y • k4d−−→ P-Y (Produits chlorés)

A• + P O•
2

k10−−→ Produits inactifs

(2.1)

(2.2)

(2.3)

(2.4)

(2.5)

(2.6)

(2.7)

(2.8)

(2.9)

(2.10)

(2.11)

(2.12)

Ce schéma, couplé à la diffusion des espèces mobiles, permet de calculer les profils de
concentration dans l’épaisseur du matériau.

Le cœur du modèle est un schéma cinétique d’autoxydation radicalaire en chaîne
(initiation, propagation, ramification, terminaison). Ce schéma, impliquant 12 réac-
tions principales, est couplé à la diffusion 1D des espèces mobiles (O2, DOC, Anti-
oxydant) via un système d’équations aux dérivées partielles (EDP). L’oxydation se ma-
nifeste par l’apparition de produits carbonylés ([CO]), comme illustré en Figure 2.4.

FIGURE 2.4 – Formation de produits carbonylés ([CO]), signature de l’oxydation de la
matrice, dans un film de PE exposé au dioxyde de chlore. Adapté de Colin et al. [2].

4

Niveau 2 : Coupures de Chaînes et Chute de la Masse Molaire

La conséquence la plus dommageable de l’oxydation à l’échelle moléculaire est la
coupure des chaînes macromoléculaires (scissions, S), principalement issue de la dé-
composition des hydroperoxydes (réactions k1u et k1b). Cette accumulation de défauts
est directement liée à la chute de la masse molaire moyenne en poids (Mw) par les
équations de Saito :

1
Mw(z, t)

−

1
Mw0

=

S(z, t)
2

− X(z, t)

(2.13)

où X est le nombre de réticulations. La chute de Mw est l’indicateur quantitatif de
l’endommagement moléculaire.

Niveau 3 : Critère de Rupture et Fragilisation Mécanique

La dégradation à l’échelle moléculaire (baisse de Mw) a une conséquence directe sur
les propriétés mécaniques à l’échelle macroscopique. En dessous d’une masse molaire
critique (Mw,crit), le réseau de chaînes enchevêtrées dans la phase amorphe perd sa
capacité à se déformer plastiquement, ce qui conduit à une fragilisation brutale du
matériau. Cette transition ductile-fragile est illustrée par la chute de l’allongement à la
rupture (εr) en fonction de Mw, comme le montre la Figure 2.5.

Le modèle couple cette fragilisation à la durée de vie via une loi de fluage empirique
qui prédit le temps à la rupture (tf ) en fonction de la contrainte, de la température, et
de la valeur de Mw à une profondeur critique.

log(tf ) = A0 +

Happ
2.3RT

+ αm log(Mw(zcrit, t)) − C log(σ)

(2.14)

Cette équation finale lie tous les niveaux : la chimie (qui génère les scissions S) contrôle
la masse molaire Mw, qui à son tour contrôle le temps à la rupture mécanique tf .

FIGURE 2.5 – Illustration de la transition ductile-fragile : chute de l’allongement à la
rupture (εr) lorsque la masse molaire (Mw) passe en dessous d’un seuil critique. Adapté
de Fayolle et al. [5].

5

2.1.3 Modèles Chemo-Mécaniques pour l’Eau de Javel

Des modèles plus récents ont tenté de simuler spécifiquement l’effet de l’eau de
Javel (HOCl). Par exemple, Tripathi et al. [7] ont développé un modèle pour la fissura-
tion sous contrainte (SCC) du PEHD. Leur approche, bien que plus axée sur la méca-
nique de la fracture, confirme la nécessité de coupler un modèle cinétique de diffusion-
réaction pour l’oxydant avec un modèle constitutif pour le comportement mécanique.
Ils montrent que la dégradation chimique accélère la propagation de fissure (Figure
2.6).

FIGURE 2.6 – Simulation de courbes contrainte-durée de vie pour du PEHD, montrant
l’accélération de la rupture en présence d’eau de Javel. Adapté de Tripathi et al. [7].

2.2 Points Durs et Piste de Recherche

2.2.1 Les Verrous Scientifiques

1. Transposition du grade de PE : Les paramètres du modèle de Colin sont spéci-
fiques au PE80. Le PE100, avec sa microstructure et ses additifs différents, a des
propriétés de transport et de réaction qui peuvent varier significativement.

2. Complexité de la chimie de l’HOCl : Contrairement au DOC, l’HOCl est en
équilibre avec l’ion ClO− et peut générer des espèces radicalaires (Cl•, HO•).
La nature exacte de l’espèce réactive qui diffuse dans le PE et sa cinétique de
réaction sont des questions ouvertes majeures.

2.2.2 Piste de Recherche et Objectifs du Stage

Face à ces verrous, la piste de ce stage est une confrontation systématique du mo-

dèle existant aux données expérimentales. Cette démarche a pour but de :
— Implémenter numériquement le modèle de Colin et al. en Python.

6

— Le valider en reproduisant les résultats de la littérature pour le PE80/DOC.
— Le confronter à des données SUEZ sur le système PE100/HOCl, en utilisant une

approche de "DOC effectif" comme première approximation.

— Quantifier les accords et les désaccords pour identifier les limites du modèle et

guider les futures améliorations nécessaires.

7

Chapitre 3

Démarche Scientifique et Outils
Numériques

Ce chapitre détaille la méthodologie employée pour simuler la dégradation du PE.
Nous présentons la formalisation mathématique du modèle, en explicitant le système
d’équations et en justifiant physiquement les conditions aux limites. Ensuite, nous dé-
crivons son implémentation numérique, en mettant l’accent sur les défis liés à la rai-
deur du système et les choix de solveurs numériques.

3.1 Formalisation Mathématique du Modèle

3.1.1 Système d’Équations de Réaction-Diffusion

Le modèle simule l’évolution spatio-temporelle de la concentration Ci(z, t) de 12
espèces chimiques. Pour chaque espèce i, son évolution est gouvernée par une équa-
tion aux dérivées partielles (EDP) qui couple la diffusion Fickienne (pour les espèces
mobiles) et les termes de réaction :

∂Ci(z, t)
∂t

= Di

∂2Ci(z, t)
∂z2

+ Ri(C(z, t))

(3.1)

où Di est le coefficient de diffusion de l’espèce i (nul si l’espèce est immobile, comme
les radicaux) et Ri(C) est un terme de source chimique non-linéaire qui représente la
somme algébrique des vitesses de toutes les réactions du schéma cinétique (Eq. 1.1 à
1.12) produisant ou consommant l’espèce i.

8

3.1.2 Description du Problème Physique et Conditions aux Limites

(CL)

La résolution de ce système d’EDP requiert la définition de conditions aux limites
(CL) qui décrivent les interactions du matériau avec son environnement. Deux cas sont
considérés : le tuyau en service et le film mince immergé en laboratoire.

Cas du Tuyau en Service

Un tuyau en service possède une interface interne (z = 0) en contact avec l’eau

désinfectée et une interface externe (z = L) en contact avec l’air.

— Interface Eau/Polymère (z = 0) :

– Oxygène et Désinfectant : On suppose que la concentration à la surface du po-
lymère est constamment à l’équilibre avec la concentration dans l’eau. C’est
une condition de Dirichlet : Ci(0, t) = Ci,eq.

– Antioxydant : Le flux d’extraction par l’eau est proportionnel à la concentra-
|z=0 = β0CAH(0, t).

tion de surface. C’est une condition de Robin : −DAH

∂CAH
∂z

— Interface Air/Polymère (z = L) :

– Oxygène : La concentration est fixée à l’équilibre avec l’air (Dirichlet).
– Désinfectant (DOC) : Le DOC est non-volatil. Son flux à travers cette interface

est donc nul. C’est une condition de Neumann : ∂CDOC

∂z

|z=L = 0.

– Antioxydant : Le flux d’évaporation est modélisé par une condition de Robin :

−DAH

∂CAH
∂z

|z=L = βLCAH(L, t).

Cas du Film Mince Immergé

Pour les essais en laboratoire, un film mince est immergé dans la solution. Les deux
faces (z = 0 et z = L) sont donc en contact avec le même environnement. Les condi-
tions sont symétriques : les CL de l’interface Eau/Polymère s’appliquent aux deux
frontières.

3.2

Implémentation et Résolution Numérique

3.2.1 Discrétisation Spatiale : la Méthode des Lignes

Le modèle a été implémenté en Python. Pour résoudre le système d’EDP, nous utili-
sons la méthode des lignes. Le domaine spatial [0, L] est discrétisé en nz points, formant
une grille d’espacement ∆z. Sur cette grille, l’opérateur de diffusion est approché par

9

une formule de différences finies centrées d’ordre 2 :

∂2Ci
∂z2

(cid:12)
(cid:12)
(cid:12)
(cid:12)
(cid:12)j

≈

Ci,j+1 − 2Ci,j + Ci,j−1
(∆z)2

(3.2)

où j est l’indice du nœud spatial. Cette démarche transforme le système d’EDP en
un grand système couplé de 12 × nz Équations Différentielles Ordinaires (EDO) en
temps, de la forme matricielle dY
dt = f (t, Y), où Y est le vecteur contenant toutes les
concentrations à tous les points de la grille.

3.2.2 Intégration Temporelle et Gestion de la Raideur

Le Défi de la Raideur Numérique

Le système d’EDO obtenu (Eq. 3.2) est numériquement "raide" (stiff). Cette pro-
priété, commune dans les problèmes de cinétique chimique, signifie que les constantes
de temps caractéristiques des différents processus (certaines réactions radicalaires très
rapides vs. la diffusion très lente) s’étalent sur plusieurs ordres de grandeur. Cela pose
un défi majeur pour l’intégration temporelle.

Limites des Schémas Explicites (ex : Runge-Kutta)

Un schéma d’intégration explicite calcule l’état du système au temps tn+1 en se

basant uniquement sur l’état et sa dérivée au temps tn.

La méthode la plus simple, celle d’Euler explicite, s’écrit : yn+1 = yn + ∆t · f (tn, yn).

Cependant, sa faible précision (erreur en O(∆t)) la rend peu pratique.

Les méthodes de Runge-Kutta (RK) améliorent considérablement la précision en
évaluant la pente en plusieurs points intermédiaires de l’intervalle de temps. Par exemple,
la méthode Runge-Kutta d’ordre 2 (RK2), ou méthode du point milieu (Figure 3.1), uti-
lise une prédiction à mi-pas pour améliorer l’estimation :

y∗ = yn +

∆t
2

f (tn, yn)

yn+1 = yn + ∆t · f

tn +

!

∆t
2

, y∗

(3.3)

(3.4)

La méthode la plus courante est celle de Runge-Kutta d’ordre 4 (RK4), qui offre
un excellent compromis entre précision (erreur en O(∆t4)) et coût de calcul pour les

10

FIGURE 3.1 – Illustration du principe d’une méthode de Runge-Kutta d’ordre 2 (RK2),
un schéma explicite performant pour les problèmes non-raides. Adapté du support de
cours MASC-FEM [6].

problèmes non-raides. Elle calcule une moyenne pondérée de quatre pentes :

k1 = f (tn, yn)

k1)

, yn +

k2 = f (tn +

∆t
2
∆t
2
k4 = f (tn + ∆t, yn + ∆tk3)

∆t
2
∆t
2

k3 = f (tn +

, yn +

k2)

yn+1 = yn +

∆t
6

(k1 + 2k2 + 2k3 + k4)

(3.5)

(3.6)

(3.7)

(3.8)

(3.9)

Cependant, pour un système raide, la stabilité des schémas explicites impose un
pas de temps ∆t extrêmement faible, dicté par la constante de temps la plus rapide du
système. Cela rendrait le calcul pour des durées de vieillissement réelles numérique-
ment non viable.

Nécessité des Schémas Implicites : BDF et Radau

Pour surmonter la raideur, il est impératif d’utiliser des schémas d’intégration ×
implicites. Ces méthodes calculent l’état yn+1 en résolvant une équation qui dépend de
yn+1 lui-même, ce qui les rend inconditionnellement stables (A-stables).

La méthode BDF (Backward Differentiation Formula) est une famille de méthodes
multi-pas implicites très efficaces pour les problèmes raides. La formule d’ordre 1,
aussi connue comme Euler implicite, s’écrit :

yn+1 = yn + ∆t · f (tn+1, yn+1)

(3.10)

11

La méthode Radau est un type de schéma de Runge-Kutta implicite, également très
performant.

Notre implémentation utilise la fonction scipy.integrate.solve_ivp, qui offre
un accès direct à ces deux solveurs robustes, en ajustant les tolérances relatives (rtol) et
absolues (atol) pour optimiser le compromis précision/temps de calcul.

3.2.3 Calcul des Grandeurs d’Intérêt

Une fois la résolution temporelle effectuée, plusieurs grandeurs d’intérêt sont cal-

culées en post-traitement à partir des profils de concentration Ci(z, t).

— Mise à jour des paramètres : Les constantes cinétiques ki et de diffusion Di
sont mises à jour pour la température de simulation T (en Kelvin) via la loi
d’Arrhenius :

k(T ) = A exp

(cid:19)

(cid:18) −Ea
RT

(3.11)

où A est le facteur pré-exponentiel, Ea l’énergie d’activation et R la constante
des gaz parfaits.

— Indice de Temps d’Oxydation (OIT) : L’OIT est supposé proportionnel à la
concentration restante en groupements antioxydants fonctionnels. Il est calculé
par la relation suivante [3] :

OIT (z, t) = ti0 ×

[AH](z, t)
[AH]0

(3.12)

où ti0 est l’OIT initial du matériau et [AH]0 sa concentration initiale. Pour les
comparaisons avec les données expérimentales macroscopiques, une valeur moyenne
d’OIT sur l’épaisseur est calculée.

3.3 Matériaux et Méthodes Expérimentales

Les simulations numériques présentées dans ce travail sont confrontées à un jeu
de données expérimentales obtenues à partir de films de polyéthylène haute densité
(PE100). Ce paragraphe détaille le matériau fourni par SUEZ et les protocoles expéri-
mentaux mis en œuvre au sein du laboratoire PIMM pour caractériser son vieillisse-
ment.

3.3.1 Matériau et Préparation des Échantillons

Le matériau étudié est un grade de PE100 fourni par SUEZ. Des films d’une épais-
seur nominale de 0,4 mm ont été préparés par compression à chaud à 180 ◦C à partir

12

de 5.6 g de granulés. Ces films constituent les échantillons de base pour les essais de
vieillissement.

3.3.2 Conditions de Vieillissement

Deux types de vieillissement accéléré ont été menés en parallèle sur une durée de 9

mois à une température constante de 40 ◦C :

— Vieillissement en eau pure : Les films sont immergés dans de l’eau déminérali-

sée.

— Vieillissement en eau de Javel : Les films sont immergés dans une solution

d’hypochlorite de sodium (HOCl).

3.3.3 Techniques de Caractérisation

Calorimétrie Différentielle à Balayage (DSC)

Les analyses DSC ont été réalisées sur un appareil TA Instruments DSC2500 à partir

d’échantillons d’environ 10 à 15 mg prélevés dans les films vieillis.

Mesure de l’Indice de Temps d’Oxydation (OIT). L’OIT, qui quantifie la stabilité
thermo-oxydative résiduelle du matériau, a été mesuré selon le protocole suivant : 1)
Rampe de chauffe de la température ambiante à 190,00 ◦C à une vitesse de 20,00 ◦C/min ×
sous atmosphère inerte (azote). 2) Palier isotherme de 5 minutes pour assurer l’équi-
libre thermique. 3) Commutation du gaz vers l’oxygène pur et maintien de l’isotherme
à 190 ◦C pendant 120 minutes. L’OIT est le temps écoulé entre l’introduction de l’oxy-
gène et le début du pic exothermique de l’oxydation.

Détermination de la Cristallinité (Xc). Le taux de cristallinité a été déterminé à partir
de l’enthalpie de fusion mesurée lors d’une rampe de chauffe de 25,00 ◦C à 220,00 ◦C
à une vitesse de 10,00 ◦C/min. Cette valeur de Xc est un paramètre d’entrée important
pour le modèle numérique.

Spectroscopie Infrarouge à Transformée de Fourier (FTIR)

Les analyses FTIR ont été menées avec un instrument PerkinElmer Frontier FT-IR
en mode de réflectance totale atténuée (ATR) avec un cristal de Diamant/ZnSe. Les
spectres ont été enregistrés sur une gamme de 4000 à 650 cm−1, avec une résolution
de 4 cm−1 et une accumulation de 24 balayages par mesure. L’objectif était de quanti-
fier l’apparition de la bande d’absorption des groupements carbonyles (C=O), située
autour de 1720-1740 cm−1.

13

Essais de Traction Mécanique

Pour évaluer l’impact du vieillissement sur les propriétés mécaniques, des essais de
traction ont été réalisés sur les films vieillis à différents temps d’exposition. Ces essais
permettent de mesurer la courbe contrainte-déformation et d’en déduire des grandeurs
clés telles que l’allongement à la rupture (εr), un indicateur particulièrement sensible à
la fragilisation du matériau.

3.4 Stratégie de Confrontation Modèle/Expérience

L’objectif final de ce travail est de confronter les prédictions du modèle cinétique
avec des données expérimentales, comme le synthétise la Figure 3.2. La démarche
adoptée est une approche par étapes, conçue pour valider ou invalider progressive-
ment les différentes briques du modèle.

FIGURE 3.2 – Approche globale du projet, montrant l’interconnexion entre les données
d’entrée (A), les blocs de simulation (B), les prédictions (C) et les données de référence
pour la validation (D).

La stratégie de confrontation s’articule comme suit :

1. Validation de l’Implémentation Numérique : La première étape, non détaillée
dans les résultats mais cruciale, a été de valider l’implémentation en Python du

14

modèle de Colin et al. en reproduisant les résultats publiés pour le système de
référence PE80/DOC [3]. Cette étape assure que notre outil de calcul est correct.

2. Découplage des Phénomènes - Perte Physique : Le vieillissement des films de
PE100 en eau pure est étudié en premier. Ce cas permet d’isoler les mécanismes
de perte physique d’antioxydant (diffusion et extraction) sans l’interférence de
la dégradation chimique. La confrontation des simulations (avec [DOC]=0) avec
les données d’OIT expérimentales permet de valider les paramètres de transport
de l’antioxydant (DAH, β0) pour le système PE100 en film mince.

3. Évaluation de la Dégradation Chimique - Approche "DOC Effectif" : Pour
le vieillissement en eau de Javel, en l’absence d’un schéma cinétique dédié à
l’HOCl, une approche exploratoire est utilisée. Une concentration "DOC effec-
tive" est introduite dans le modèle. Cette approche permet une première évalua-
tion qualitative de la capacité du modèle à simuler une dégradation accélérée.
La confrontation des prédictions (OIT, [CO]) avec les données expérimentales
SUEZ permet d’analyser les accords et, surtout, les désaccords, qui sont riches
d’enseignements sur les différences de réactivité entre DOC et HOCl.

4. Validation du Couplage Mécano-Chimique : La dernière étape consiste à corré-
ler les indicateurs de dégradation chimique (chute de l’OIT) avec les indicateurs
de dégradation mécanique (chute de l’allongement à la rupture). Cette confron-
tation vise à vérifier expérimentalement la chaîne de causalité postulée par le
modèle, où la dégradation chimique gouverne la fragilisation et donc la durée
de vie du matériau.

Cette démarche progressive permet de tirer des conclusions robustes à chaque étape
et d’identifier précisément les forces et les faiblesses du modèle lorsqu’il est appliqué
en dehors de son domaine de validation initial. Les résultats de cette confrontation
sont présentés dans le chapitre suivant.

15

Chapitre 4

Résultats et Analyses

Ce chapitre présente les résultats obtenus par l’implémentation numérique du mo-
dèle cinétique. La première section détaille la stratégie et les paramètres de simulation.
La seconde présente une analyse des résultats préliminaires de simulation pour un cas
de référence. La troisième section est consacrée à la validation du modèle par confron-
tation directe avec les données de la littérature. Enfin, la dernière section aborde la
confrontation avec les données expérimentales SUEZ.

4.1 Stratégie de Modélisation et Paramètres de Simula-

tion

Toutes les simulations sont réalisées avec le code Python développé durant ce stage.
Sauf mention contraire, les paramètres cinétiques (constantes de vitesse ki, énergies
d’activation Ea) et de transport (Di, βi) sont ceux publiés par Colin et al. [2, 3] pour le
système de référence PE80/DOC.

4.1.1 Paramètres du Système de Référence (Validation)

Pour la validation du modèle par rapport à la littérature, les paramètres suivants

sont utilisés pour simuler un tuyau de PE80 :
— Géométrie : Épaisseur L = 4,5 mm.
— Propriétés Initiales : Masse molaire Mw0 = 150 kg/mol, concentration en sites
oxydables [PH]0 = 60 mol/L, OIT initial de référence ti0 = 165 min, et cristallinité
de référence Xc = 0.45. Cette dernière valeur est cohérente avec les mesures DSC
réalisées dans le cadre de ce stage.

— Profil d’Antioxydant Initial : Pour tenir compte de la perte de stabilisant lors
de la mise en œuvre, un profil initial non-uniforme pour [AH](z, t = 0) est im-
plémenté, reproduisant une déplétion de 25% en surface externe (OSL) [4].

16

— Conditions Opératoires : Température et concentrations en DOC spécifiées pour

chaque cas de validation.

4.2 Analyse des Profils de Concentration Simulés (Résul-

tats Préliminaires)

Une première simulation pour un cas typique (tuyau PE80 exposé à 1 ppm de DOC
à 40◦C) permet d’observer l’évolution des profils de concentration des espèces mobiles
(O2, DOC, AH) à travers l’épaisseur sur une durée de 15 ans (Figure 4.1).

FIGURE 4.1 – Profils de concentration simulés pour l’Oxygène, l’Antioxydant (AH) et
le Dioxyde de Chlore (DOC) à différents temps de vieillissement.

L’analyse de ces profils révèle plusieurs phénomènes physiques et chimiques clés

prédits par le modèle :

— Profil de [DOC] : Le désinfectant pénètre depuis la surface interne (z = 0) et est
rapidement consommé. Sa concentration devient négligeable au-delà de 0,5 mm,
indiquant que sa pénétration est limitée par sa consommation réactive.

— Profil de [AH] : L’antioxydant est consommé en surface interne par réaction
avec le DOC, créant un "front de déplétion" qui progresse dans l’épaisseur avec
le temps. Dans le cœur du matériau, sa concentration diminue plus lentement
par diffusion vers les deux interfaces.

— Phénomène Asymétrique de [O2] : Le profil d’oxygène présente un comporte-
ment non trivial. Aux temps courts (ex : 3.67 ans), un "puits" de concentration
apparaît près de la surface interne. Cela s’explique par le fait que l’initiation
de l’oxydation par le DOC consomme localement de l’oxygène plus rapidement
qu’il n’est fourni par la diffusion depuis l’eau. Ce phénomène, signature du cou-
plage réaction-diffusion, disparaît aux temps longs lorsque la consommation
d’oxygène devient limitée par la diffusion depuis les deux interfaces.

17

4.3 Validation du Modèle Numérique par Confrontation

à la Littérature

Cette étape vise à valider notre implémentation en comparant ses prédictions à des

résultats publiés.

4.3.1 Validation des Profils d’OIT (Réf. Figure 3, Colin et al. 2009)

La Figure 4.2 confronte les profils d’OIT simulés (gauche) avec la figure de référence
de la littérature (droite). L’accord qualitatif est très bon : les plateaux d’OIT, les fronts de
déplétion et leur évolution avec la concentration en DOC sont bien reproduits. L’écart
quantitatif sur la position exacte des fronts est attribuable aux incertitudes sur les pa-
ramètres cinétiques et de transport du modèle de base.

(a) Simulation

(b) Littérature [3]

FIGURE 4.2 – Confrontation des profils d’OIT simulés (a) et expérimentaux (b) pour 5
scénarios de vieillissement à 40◦C.

4.3.2 Validation des Profils de Carbonyles et d’OIT (Réf. Figure 4,

Colin et al. 2009)

La Figure 4.3 confronte les profils simulés d’OIT et de [CO] avec ceux de la littéra-
ture pour un cas de vieillissement sévère. L’accord est excellent et valide le couplage
clé du modèle : la formation de produits d’oxydation ([CO]) est correctement prédite
comme étant confinée à la zone de surface où l’antioxydant (OIT) a été totalement
consommé.

Cette double validation confirme que notre implémentation numérique est fonc-
tionnelle et apte à être utilisée pour la confrontation avec les nouvelles données expé-
rimentales.

18

(a) Simulation

(b) Littérature [3]

FIGURE 4.3 – Confrontation des profils d’OIT (pointillés) et de Carbonyles (symboles
pleins) simulés (a) et expérimentaux (b) après 99 jours à 40◦C.

4.4 Confrontation avec les Données SUEZ (Films PE100)

4.4.1 Vieillissement en Eau Pure (H2O) : Cinétique de Perte Physique

L’étude du vieillissement en eau pure permet d’isoler la perte physique d’antioxy-
dant. L’OIT simulé, présenté en Figure 4.4, est calculé comme la moyenne du profil
‘OIT(z,t)‘ sur l’épaisseur du film. La Figure 4.4 montre un accord qualitatif satisfaisant.

FIGURE 4.4 – Évolution de l’OIT vs. Temps pour films PE100 en eau pure à 40 ◦C.
Confrontation : Données SUEZ (bleu) et Simulation (rouge).

Quantitativement, le modèle prédit une perte d’OIT légèrement plus rapide que celle
observée expérimentalement. À 9 mois, la perte simulée est de 52 min contre une perte

19

expérimentale de 41,8 min. Cela suggère que les paramètres de transport (DAH, β0) is-
sus de la littérature pour des tuyaux PE80 surestiment légèrement la cinétique de perte
pour ces films minces de PE100.

4.4.2 Vieillissement en Eau de Javel (HOCl) : Dégradation Chimique

Accélérée

La confrontation en milieu HOCl (Figure 4.5) évalue l’impact chimique du désin-

fectant. L’analyse de la Figure 4.5 révèle deux points majeurs :

FIGURE 4.5 – Évolution de l’OIT (axe gauche) et de l’Indice Carbonyle (axe droit) en
fonction du temps pour des films de PE100 sous HOCl à 40 ◦C. Données SUEZ : OIT
(vert), Indice Carbonyle (orange). Simulation : OIT (violet), [CO] en surface z=0 (mar-
ron).

1. Consommation de l’antioxydant (OIT) : L’approche du “DOC effectif” repro-
duit bien la chute rapide initiale de l’OIT, mais surestime la vitesse de dégra-
dation aux temps longs, prédisant une déplétion totale plus précoce que dans
l’expérience.

2. Formation des produits d’oxydation ([CO]) : Le modèle prédit une augmenta-
tion monotone de la concentration en carbonyles en surface ([CO](z = 0, t)), tan-
dis que les données expérimentales montrent un comportement non-monotone.
Ce dernier pourrait s’expliquer par une compétition entre la formation de pro-
duits oxydés et leur lixiviation hors du film, un phénomène non inclus dans le
modèle actuel.

20

Ces résultats confirment que l’HOCl accélère significativement la dégradation. Cepen-
dant, les écarts quantitatifs soulignent les limites de l’analogie avec le DOC et la néces-
sité d’un modèle cinétique spécifique à l’HOCl pour des prédictions fines.

Analyse Spectroscopique de l’Oxydation par FTIR

L’oxydation de la matrice polyéthylène se manifeste par l’apparition de groupe-
ments fonctionnels oxygénés, notamment des cétones et des acides carboxyliques. Ces
groupements présentent une bande d’absorption caractéristique dans la région de 1700
à 1740 cm−1 due à la vibration d’élongation de la double liaison carbonyle (C=O).

La Figure 4.6 présente la superposition des spectres FTIR, obtenus en mode ATR et
normalisés sur une bande de référence, pour les films de PE100 vieillis en eau de Javel
(HOCl) à différents temps d’exposition.

FIGURE 4.6 – Évolution de la bande d’absorption carbonyle ( 1720 cm−1) dans les
spectres FTIR des films de PE100 en fonction du temps de vieillissement sous HOCl.

L’analyse de ces spectres confirme une oxydation progressive de la surface du poly-
mère. On observe l’émergence et la croissance d’une bande d’absorption large centrée
autour de 1720 cm−1. L’intensité de ce pic, qui est absente sur le matériau non-vieilli,
augmente avec la durée d’exposition, ce qui constitue une preuve qualitative directe
de l’accumulation de produits d’oxydation.

La hauteur de ce pic est donc utilisée dans la suite de ce travail comme un indica-
teur quantitatif, appelé "Indice Carbonyle", pour suivre la cinétique de la dégradation
chimique de la matrice polymère.

21

4.4.3 Impact sur les Propriétés Mécaniques : Vers une Prédiction de

Durée de Vie

Après avoir caractérisé et modélisé la dégradation chimique du PE100, l’étape finale
consiste à évaluer l’impact de cette dégradation sur les propriétés mécaniques, qui
gouvernent la durée de vie en service du matériau. L’allongement à la rupture (εr) est
un indicateur particulièrement sensible à la fragilisation induite par les coupures de
chaînes.

La Figure 4.7 met en corrélation directe l’évolution de cet indicateur mécanique

avec l’indicateur chimique (OIT) pour les deux milieux de vieillissement.

FIGURE 4.7 – Corrélation entre la dégradation chimique (OIT, axes de droite) et la dé-
gradation mécanique (Allongement à la rupture, axes de gauche) pour des films de
PE100 vieillis à 40◦C en eau pure et en eau de Javel (HOCl).

L’analyse de ce graphique révèle une corrélation très nette :
— En eau pure (H2O) : La perte d’OIT est très limitée sur 9 mois. Parallèlement,
l’allongement à la rupture reste élevé et stable, aux alentours de 950 %. L’absence
de dégradation chimique significative se traduit par une absence de fragilisation
mécanique.

— En eau de Javel (HOCl) : Le comportement est radicalement différent. On ob-
serve une chute simultanée et drastique de l’OIT et de l’allongement à la rup-
ture. Après seulement 6 mois, alors que l’OIT a chuté de près de 93 % (de 271 min
à 19 min), l’allongement à la rupture s’est effondré de plus de 80 % (passant de
≈ 930 % à moins de 200 %).

22

Ces résultats expérimentaux démontrent de manière claire la chaîne de causalité
postulée par le modèle de Colin et al. : la consommation de l’antioxydant (chute de
l’OIT) permet l’oxydation de la matrice polymère, qui induit des coupures de chaînes
(baisse de Mw), lesquelles provoquent à leur tour une fragilisation sévère du matériau
(chute de εr).

Perspective : Vers la Prédiction de Durée de Vie. Le Niveau 3 du modèle cinétique
formalise ce lien via une loi de durée de vie (Eq. 2.14) qui relie le temps à la rupture
(tf ) à la masse molaire (Mw). Les données présentées ici valident qualitativement l’hy-
pothèse fondamentale de ce couplage mécano-chimique.

La validation quantitative complète, qui consisterait à simuler les temps de rupture
pour différentes contraintes et à les comparer à des données expérimentales de rupture
sous charge, constitue une perspective de recherche directe et prometteuse, s’appuyant
sur les fondations établies par ce travail.

23

Chapitre 5

Conclusions et Perspectives

Ce stage avait pour objectif d’évaluer la pertinence d’un modèle cinétique de dé-
gradation existant, développé pour le système PE80/DOC, dans le contexte industriel
d’un grade PE100 exposé à l’hypochlorite de sodium (HOCl). Pour ce faire, le modèle
a été implémenté numériquement en Python, validé par rapport aux données de la lit-
térature, puis confronté à des données expérimentales de vieillissement fournies par
SUEZ.

Les principales conclusions de ce travail sont les suivantes :
— L’implémentation numérique du modèle s’est avérée fonctionnelle et capable de
reproduire avec une bonne fidélité les résultats de référence publiés, validant
ainsi notre outil de simulation.

— La confrontation avec les données de vieillissement en eau pure a permis de
confirmer que le modèle, avec ses paramètres de transport d’origine, capture
correctement la cinétique de perte physique d’antioxydant, bien qu’avec une
légère surestimation de la vitesse.

— Pour le vieillissement en eau de Javel, une approche exploratoire utilisant une
concentration "DOC effective" a montré que le modèle pouvait simuler qualita-
tivement l’accélération de la dégradation. Cependant, les écarts quantitatifs et
le comportement non-monotone de l’indice carbonyle expérimental soulignent
les limites de cette analogie et la nécessité de développer un schéma cinétique
spécifique à l’HOCl.

— L’analyse des propriétés mécaniques a mis en évidence une corrélation expé-
rimentale directe et très nette entre la dégradation chimique (chute de l’OIT)
et la fragilisation du matériau (chute de l’allongement à la rupture). Ce résultat
clé valide l’hypothèse fondamentale du couplage mécano-chimique sur laquelle
repose le modèle.

Perspectives. Ce travail a permis de poser des fondations solides et ouvre plusieurs
pistes de recherche prometteuses. À court terme, une étude de sensibilité sur les para-

24

mètres les plus influents (DDOC, k8d) pourrait permettre d’améliorer l’ajustement quan-
titatif du modèle. À plus long terme, le développement d’un schéma cinétique dédié à
la chimie complexe de l’HOCl est indispensable pour atteindre une capacité de prédic-
tion fine. La validation finale du modèle de durée de vie, en confrontant les temps de
rupture simulés à des essais expérimentaux sous charge, constitue l’étape ultime de ce
projet de recherche.

25

Bibliographie

[1] Xavier Colin, Ludmila Audouin, and Jacques Verdu. Towards a non empirical ki-
netic model for the lifetime prediction of polyethylene pipes transporting drinking
water. 286(1) :81–88, 2009.

[2] Xavier Colin, Ludmila Audouin, Jacques Verdu, Magali Rozental-Evesque, Benja-
min Rabaud, Florencio Martin, and Francis Bourgine. Aging of polyethylene pipes
transporting drinking water disinfected by chlorine dioxide. i. chemical aspects.
Polymer Engineering & Science, 49(7) :1429–1437, 2009.

[3] Xavier Colin, Ludmila Audouin, Jacques Verdu, Magali Rozental-Evesque, Benja-
min Rabaud, Florencio Martin, and Francis Bourgine. Aging of polyethylene pipes
transporting drinking water disinfected by chlorine dioxide. part ii—lifetime pre-
diction. Polymer Engineering & Science, 49(8) :1642–1652, 2009.

[4] Xavier Colin, Jacques Verdu, and Benjamin Rabaud. Stabilizer thickness profiles
in polyethylene pipes transporting drinking water disinfected by bleach. Polymer
Engineering & Science, 51(8) :1541–1549, 2011.

[5] B. Fayolle, E. Richaud, X. Colin, and J. Verdu. Mechanism of degradation induced
embrittlement in polyethylene. Polymer Degradation and Stability, 92(2) :231–238,
2007.

[6] Juan-Pablo Marquez Costa and Jan Neggers. Material and structural computation
by fem (masc-fem), session 7 - implicit schemes for integration of material beha-
viour laws. Support de cours, Master MAGIS, Arts et Métiers, 2023.

[7] Anshul Tripathi, Susan C. Mantell, and Jia-Liang Le. Chemo-mechanical modeling
of static fatigue of high density polyethylene in bleach solution. International Journal
of Solids and Structures, 217-218 :90–105, 2021.

[8] W. Yu, B. Azhdar, D. Andersson, T. Reitberger, J. Hassinen, T. Hjertberg, and U. W.
Gedde. Deterioration of polyethylene pipes exposed to water containing chlorine
dioxide. Polymer Degradation and Stability, 96(5) :790–797, 2011.

26


