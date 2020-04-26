function constantes=retconstantes(id);
% Constantes physiques en unit�s MKSA
% http://www.h-k.fr/publications/data/adc.ps__annexes.physique-chimie.pdf
% constantes=retconstantes :est une structure � champs de champs
%    'c'  'e'  'NA' 'G' 'R' 'F' 'kB' 'h' 'Pba' 'me' 'mn' 'mp' 'mu0' 'h_bar' 'eps0' 'Z0' 'A' 'prefixes'
% constantes=retconstantes(nom) valeur de constante.nom  (nom est une chaine de caract�res)
% 
% Vitesse de la lumi�re c=2.99792458e8 m.s-1 (la constante lumineuse demeurera desormais l� dans votre cervelle)
% Charge �l�mentaire e=1.60219e-19 C
% Nombre d�Avogadro NA=6.02204e23 mol-1
% Constante gravitationnelle G=6.672e-11 N.m2.kg-2
% Constante des gaz parfaits R=8.3144 J.K-1.mol-1
% Constante de Faraday F=96484 C.mol-1
% Constante de Boltzmann kB = 1.38066e-23 J.K-1
% Constante de Planck h = 6.62617e-34 J.s
%                     h_bar =h/(2*pi)
% Constante de Pba: Pba=hc/kB
% Masse de l��lectron me = 9.10953e-31 kg
% Masse du neutron mn = 1.675e-27 kg
% Masse du proton mp = 1.673e-27 kg
% Permittivit� du vide eps0 = 8.85419e-12 F.m-1
% Perm�abilit� du vide mu0 = 4.e-7*pi H.m-1
% constantes=retconstantes;
% ou c=retconstantes('c');
%
%%%  Exemple
% c=retconstantes('c')
%
% See also: RETDB

% constantes=struct('c',2.99792458e8,...
%   'e',1.60219e-19,...
%   'NA',6.02204e23,...
% 'G',6.672e-11,...
% 'R',8.3144,...
% 'F',96484,...
% 'kB',1.38066e-23,...
% 'h',6.62617e-34,...
% 'me',9.10953e-31,...
% 'mn',1.675e-27,...
% 'mp',1.673e-27,...
% 'mu0',4.e-7*pi,...
% 'A',1.e-10);% modif 6 2014 avec PBA
constantes=struct('c',2.99792458e8,...
'e',1.60217656535e-19,...
'NA',6.02204e23,...
'G',6.672e-11,...
'R',8.3144,...
'F',96484,...
'kB',1.380648813e-23,...
'h',6.62606957e-34,...
'me',9.10953e-31,...
'mn',1.675e-27,...
'mp',1.673e-27,...
'mu0',4.e-7*pi,...
'A',1.e-10);
constantes.Pba=constantes.h*constantes.c/constantes.kB;
constantes.h_bar=constantes.h/(2*pi);
constantes.eps0=1/(constantes.mu0*constantes.c^2);
constantes.Z0=sqrt(constantes.mu0/constantes.eps0);
constantes.prefixes=struct('yotta',1.e24,'zetta',1.e21,'exa',1.e18,'peta',1.e15,'tera',1.e12,'giga',1.e9,'mega',1.e6,'kilo',1.e3,'hecto',1.e2,'deca',1.e1,...
'deci',1.e-1,'centi',1.e-2,'milli',1.e-3,'micro',1.e-6,'nano',1.e-9,'pico',1.e-12,'femto',1.e-15,'atto',1.e-18,'zepto',1.e-21,'yocto',1.e-24);
if nargin>0;constantes=getfield(constantes,id);end;