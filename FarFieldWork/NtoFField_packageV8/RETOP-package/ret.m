%---------------------------------------------------------
% LISTE DES PROGRAMMES UTILISES COURAMMENT PAR RETICOLO 
%---------------------------------------------------------
%
%  Programmes permettant de traiter un probl�me de r�seau 
%---------------------------------------------------------
%  retep          :construction de epsilon et mu a partir de l'indice (1D)
%  ret2ep         :construction de epsilon et mu a partir de l'indice (2D)
%  retu           :construction d'un maillage de texture
%  retdisc        :discontinuit�s d'un maillage de texture
%  retcouche      :calcul du descripteur d' une texture (diagonalisation)
%  retb           :conditions aux limites en haut et en bas
%  retc           :matrice s sur une hauteur h
%  retreseau      :mise en forme des efficacit�s et des ordres dans le cas d'un reseau (thales)
%  retbragg       :reseaux holographiques en incidence conique
%
% 
%   Mie, Polarisabilit�, YLM
% ---------------------------------------------------------
%  retMie         : theorie des spheres de Mie
%  ret_bessel_spherique: fonctions de bessels spheriques
%  retlegendre    : fonctions de legendre
%  retexclucion   : calcul de la matrice L associ�e � un domaine d'exclusion
%  retcarminati   : particules dans un systeme stratifie eventuellement periodique
%  ret_plusieurs_sources: Flux total diffract� par plusieurs sources
%  ret_Polarisabilite_cylindrique: calcul de la matrice polarisabilit� d'un objet de revolution d'axe oz
%  ret_rotation_polarisabilite: transformation d'un multipole et de la matrice de polarisation associ�e dans une rotation d'espace
%  ret_theoreme_addition_YLM: theoreme d'addition des YLM
%  ret_inc2YLM    : decomposition d'un champ incident produit par une onde plane ou une source sur les harmoniques de Mie r�gulieres
%  retcalP_bonod  : Calcul de la matrice de polarisabilite 'Bonod' d'un systeme spherique ou cylindrique
%  retYLM: calcul des harmoniques (multipoles)
%  ret_bonod2multipole:construction des multisources � partir d'une matrice de polarisabilit� bonod
% retellipsoide:   ellipsoide electrostatique
%
%  Utilitaires de manipulation des matrices S (ou G)
% ---------------------------------------------------------
%  retgs          :transformation S a G ou T dans tous les sens
%  retrenverse    :renversement matrice S 
%  rettr          :translation matrice S (translation objet)
%  retstr         :matrice S d'une translation 
%  rets1          :matrice S unit�
%  retss          :produit matrices S (ou G)
%  retsp          :puissance matrices S (ou G)
%  rettronc       :troncature des matrices S
%  reteval        :�valuation d'�l�ments d'une matrice S
%  retauto        :calcul des matrices S et des tab associ�es
%  retrands       :generation de matrices 2 2 n'amplifiant pas l'energie
%  retpassage     :passage entre deux milieux avec des PML reelles differentes.
%  retrenverse    :renversement matrice S 
% 
%   Divers
% ---------------------------------------------------------
%  retchamp       :calcul des champs E et H
%  rettronc       :troncature des matrices S
%  retfp          :Fabry Perot generalise
%  retfps         :recherche des r�sonances dans une couche invariante (source interne)
%  retfpm         :recherche des r�sonances dans une couche invariante (mode interne)
%  reti           :couche inclinee(1D)
%  rets           :sources ponctuelles
%  retinc         :champ incident consid�re comme source
%  retds          :diagonalisation matrice S
%  retbloch       :modes de Bloch 
%  rettbloch      :diagrammes de dispersion avec modes de Bloch 
%  retbrillouin   :modes de Brillouin
%  retmaystre     :diagrammes de dispersion en cherchant les p�les
%  retvg          :calcul de vg et GVD � partir de points pris sur des diagragmes de dispersion
%  retop          :d�veloppement en ondes planes
%  retcaoo        :changement de coordonn�es
%  retneff        :indice effectif d'une couche
%  retautomatique :g�n�ration automatique de tranches a partir d'un u 2D, instalation de PML
%  rethomogenise  :homog�n�isation d'un maillage de texture ou de tron�on
%  retrotation    :rotation d'un maillage de texture de pi/2 -pi/2  ou pi
%  redecoupe      :decoupe dans un maillage de tron�on de la partie comprise entre 2 hauteurs
%  retpoynting    :calcul du vecteur de Poynting
%  retvm          :calcul du volume modal
%  retpoint       :champ d'une source ponctuelle dans un milieu homog�ne isotrope
%  retfaisceau    :onde plane faisceau gaussien
%  retoc          :champ d'une source ponctuelle sur un dioprte plan ( onde cylindrique )
%  retsommerfeld  :demi plan de Sommerfeld
%  retjurek       :analyse des singularites aux aretes
%  retq           :calcul de q a partir des p�les en h
%  retmode        :modes d'un empilement de couches ou d'un objet 2 D
%  retpmode       :poynting et lorentz d'un mode symetrique(fente dans du metal ..)
%  retmarcuse     :modes d'une fibre circulaire
%  retabeles      :R et T d'un empilement (possibilit� de variable vectotielles)
%  retyeh         :empilement de couches anisotropes
%  retfr          :calcul de diffraction de Fresnel
%  retbahar       :�mission d'une source 1 D  ou 2 D dans un milieu stratifie (0 D ou cylindrique). Pour les plans calcul du champ et DOP
%  retindice      :indice de materiaux (or ag ...) en fonction de lambda et aussi mode fondamental d'une fibre circulaire
%  retdrude       :fit d'un modele de Drude � partir de retindice
%  retplasmon     :constante de propagation d'un plasmon 
%  retpl          :recherche du plasmon dans un champ par produit scalaire
%  ret_int_lorentz:integrale intervenant dans la relation de lorentz
%  retrodier      :interpolation de r et t au bord de gap
%  retdb          :conversion en db/mm d'un indice complexe
%  retconstantes  :Constantes physiques
%  retw1          :construction des maillage de textures associes au W1
%  retfpcomplexe  :calcul des min et max de (abs(X-X1)*abs(X-X2))^2
%  ret_divise_cao :partition du changement de coordonnees ( pour objets symetrises)
%  rethelp_popov  :help des objets de revolution methode Popov
%  retcr          :help des objets de revolution methode cylindrique radiale.Transformation en coordonn�es Popov
%  retcr_centre   :source au centre (methode cylindrique radiale)
%  retcr_op       :diagramme de rayonnement d'une source exentr�e (methode cylindrique radiale)
%  retdisque      :etude de certains objets de r�volution par la methode cylindrique radiale(antenne de Benjamin ,Isabelle ..)  
%  retbolomey     :equation integrale volumique de Bolomey
%  retmondher     :transformation d'un champ reticolo en MKSA
%  retgraphene    : conducticite du graphese, matrices S associees
%  ret_modes_Benjamin:  Amplitude �mise sur un mode par des sources dansun systeme stratifi�
%  ret_Lindhar_Mermin: formules de Lindhar Mermin  fonctions fl et ft de GW Ford  WH Weber
%   ret_micro_structured_fiber : fibres micro structur�es
% 
%  Test et v�rifications
% ---------------------------------------------------------
%  retener        :test de conservation de l'�nergie
%  rettestobjet   :v�rification de l'objet construit ( mais pas uniquement utilis� pour les tests: permet le calcul de epsilon mu de l'objet) 
%  rettmode       :trace du champ d'un mode d'une couche
%  retmaxwell     :v�rification �quation de maxwell
% 
%  Utilitaires graphiques
% ---------------------------------------------------------
%  retcolor       :idem a pcolor mais complet...
%  retplot        :trac� de courbes en cours de calcul
%  retfleche      :trac� de fl�ches
%  retsubplot     :subplot avec option 'align' compatible matlab 6
%  retget         :r�cup�ration de points sur une courbe, �limination de points,coupes d'une image
%  rettetramesh   :comme tetramesh mais vue '�clat�e'
%  ret_vol_maillage:Mise en ordre d'un maillage et calcul eventuel de volumes et surfaces
%  retwarp        :plaquage d'une image sur une surface
%  ret_normalise_DOP  :renormalisation d'un DOP 3D noir
% 
%  Retmatlab
% ---------------------------------------------------------
%  rettf          :tf d'une fonction avec FFT
%  retversion     :num�ro de la version Matlab
%  retsqrt        :sqrt avec d�termination des C.O.S
%  retsparse      :transforme les tableaux d'un cell-array en sparse
%  reteig         :diagonalisation m�me pour matrices sparses
%  retinterp      :interpolation  avec points �gaux
%  retinterps     :interpolation des matrices S 2x2
%  retfind        :find a et ~a
%  retcontour     :rafinement par cadillac de contour
%  retmax         :max d'un tableau a plusieurs dimensions
%  retprod        :produit de matrices dont certaines sont du type {real,imag}
%  retsave        :stockage securise sur deux fichiers
%  retfig         :ouverture et stockage des figures
%  retefface      :efface fichiers temporaires ou figures
%  retprofiler    :camembert du profiler
%  retbessel      :Acc�l�ration du calcul des fonctions de bessel d'indice entier<0
%  retwhos        :construction d'une structure contenant les variables
%  retind2sub     :identique � ind2sub de matlab mais plus rapide ( facteur 1.5)
%  retsub2ind     :identique � sub2ind de matlab mais plus rapide ( facteur 3)
%  retlog1p       :identique � log1p mais plus rapide et marche pour les complexes
%  retexpm1       :tr�s preferable � expm1
%  retslash       : / apres conditionnement
%  retantislash   :\ apres conditionnement
%  retsize        : size d'un ensemble d'objets
%  retcaxisequal  : egalise les caxis de tous les subplots de 2 figures
%  retreshape     : reshape sans erreur avec matrices creuses
%  retatan2       : atan2 pour des complexes
% 
%  Utilitaires g�n�raux
% ---------------------------------------------------------
%  retcauchy      :d�tection des poles et coupures d'une fonction complexe
%  retcadilhac    :z�ros fonction complexe
%  retzp          :tous les z�ros et p�les d'une fonction complexe dans un domaine
%  retrecuit      :recuit simul�
%  retgauss       :int�gration m�thode de gauss
%  retcombine     : mixture des discontinuites de retgauss
%  retfresnel     :fonction de Fresnel
%  retf           :calcul des cooefficients de Fourier
%  rettfpoly      :tf de polygones 
%  rettfsimplex   :tf de triangles et t�tra�dres (et d�riv�es) 
%  retintegre     :int�gration de y(x)*exp(g*x) quand y(x) a des p�les voisins de l'axe reel
%  retderivee     :d�rivation d'une fonction donnee en des points
%  retheure       :heure  minutes secondes
%  retcompare     :diff�rence des valeurs num�riques de 2 cell-array
%  retlorentz     :�chantillonnage d'une lorentzienne
%  retextremum    :recherche des extremum d'une fonction definie par des points x,y
%  retsum         :c=a+b m�me si a=[] ou b=[]
%  retsort        :comme sort mais donne aussi l'ordre inverse
%  retrand        :engendre n nombres qui tir�s au hasard donnent une loi donnee
%  rethist        :calcul de la densite de probabilit�
%  retscale       :c=sum(a.*b...  .*d) 
%  rettexte       :ecriture de donnees numeriques (structures sous forme 'developpee')
%  retelimine     :elimination  des elements egaux 
%  retassocie     :association 2 a 2, de 2 ensembles de vecteurs tres voisins 2 a 2 mais dans le desordre
%  retpermute     :permutation d'elements
%  retcolonne     :transforme un tableau en colonne ou ligne
%  retconvS       :conversion des matrices S 2 2
%  retfullsparse  :transformation en matrice sparse ou full suivant la densite
%  retdiag        :spdiags(a,0,n,n)  ..etc..
%  retapod        :apodisation
%  retio          :creation et utilisation de fichiers temporaires ou permanents 
%  retv6          :transformation de fichiers en version 6
%  retcal         :calcul non deja fait ?
%  retx           :generation d'une nouvelle valeur de xx 
%  retdeplacement :rotations translation en 2D et 3D 
%  retbornes      :limitation a une boite
%  retprecis      :calculs en precision etendue
%  retpoids       :poids periodiques pour les elements finis
%  rettic,rettoc  :comme tic,toc mais donne le temps cpu 
%  retfont        :pour une font de qualite aussi traces multiples avec divers types de lignes
%  retsphere      :maillage sphere de rayon 1 pour integration en surface
%  retcardan     : formule de cardan vectorialis�e
%  ret            :clear,efface fichiers temporaires et figures,option 3 a retss,coupure de Maystre
%  retexemple     :gabaris pour �crire des programmes
%
%  Per il divertimento:
%---------------------------------------------------------
% retsudoku      : resolution des sudoku de toutes dimensions
% retsamourai    : resolution des samourai
%
% See also: 
%  RET            RETABELES               RETATAN2       RET_BESSEL_SPHERIQUE    RET_DIVISE_CAO     RET_INT_LORENTZ     RET_VOL_MAILLAGE   RET2EP      
%  RETANTISLASH   RETAPOD                 RETASSOCIE         RETAUTO             RETAUTOMATIQUE     RETB           
%  RETBAHAR       RETBESSEL               RETBLOCH           RETBOLOMEY          RETBORNES          RETBRAGG           
%  RETBRILLOUIN   RETC                    RETCAUCHY          RETCADILHAC         RETCAL             RETCALP_BONOD      RETCARDAN   RETCAOO
%  RETCARMINATI   RETCAXISEQUAL           RETCHAMP           RETCOLONNE          RETCOMBINE         RETCOLOR
%  RETCOMPARE     RETCONSTANTES           RETCONTOUR         RETCONVS            RETCOUCHE          RETCR              RETCR_CENTRE
%  RETCR_OP       RETDB                   RETDECOUPE         RETDEPLACEMENT      RETDERIVEE         RETDIAG
%  RETDISC        RETDISQUE               RETDRUDE           RETDS               RETEFFACE          RETEIG
%  RETELIMINE     RETELLIPSOIDE           RETENER            RETEP               RETEVAL            RETEXEMPLE         RETEXCLUSION    
%  RETEXPM1       RETEXTREMUM             RETF               RETFAISCEAU         RETFIG             RETFIND        
%  RETFLECHE      RETFONT                 RETFP              RETFPCOMPLEXE       RETFPM             RETFPS         
%  RETFR          RETFRESNEL              RETFULLSPARSE      RETGAUSS            RETGET             RETGRAPHENE
%  RETGS          RETHELP_POPOV           RETHEURE           RETHIST             RETHOMOGENISE      RET_LINDHAR_MERMIN   RET_THEOREME_ADDITION_YLM   RETI
%  RETINC         RET_INC2YLM             RETIND2SUB         RETINDICE           RETINTEGRE         RETINTERP           RETINTERPS
%   RETIO          RETJUREK                RETLEGENDRE        RETLOG1P            RETLORENTZ         RETMARCUSE     RET_MICRO_STRUCTURED_FIBER
%  RETMAX         RETMAXWELL              RETMAYSTRE         RETMIE              RETMODE            RET_MODES_BENJAMIN  RETMONDHER
%  RETNEFF        RET_NORMALISE_DOP       RETOC              RETOP               RETPASSAGE         RETPERMUTE         RETPL
%  RET_PLUSIEURS_SOURCES                  RETPMODE           RETPOIDS            RET_POLARISABILITE_CYLINDRIQUE
%  RETPOINT       RETPOYNTING             RETPRECIS          RETPROD             RETPROFILER        RETQ
%  RETRAND        RETRANDS                RETRECUIT          RETRENVERSE         RETRESEAU          RETRODIER
%  RETROTATION    RET_ROTATION_POLARISABILITE                RETS                RETS1              RETSPHERE           RETSAVE
%  RETSCALE       RETSIZE                 RETSLASH           RETSOMMERFELD       RETSORT            RETSP
%  RETSPARSE      RETSQRT                 RETSS              RETSTR              RETSUB2IND         RETSUBPLOT
%  RETSUM         RETTBLOCH               RETTESTOBJET       RETTETRAMESH        RETTEXTE           RETTF
%  RETTFPOLY      RETTFSIMPLEX            RETTMODE           RETTR               RETTRONC           RETU 
%  RETV6          RETVERSION              RETVG              RETVM               RETW1              RETWARP
%  RETWHOS        RETX                    RETYEH             RETZP          
% 

clear;retio;retss(3);retfig;colordef white;clear retsqrt retep ret2ep;format compact;
set(0,'DefaultAxesLineStyleOrder','default','DefaultAxesColorOrder',[[0,0,0];[1,0,0];[0,1,0];[0,0,1];[0,0.5,0.5];[0.5,0.5,0];[0.5,0,0.5]]);
sparms=spparms;if length(sparms)>8;sparms(11)=0.01;spparms(sparms);end;% parametre pour la division des matrices creuses
rand('state',sum(100*clock));% etat de rand




