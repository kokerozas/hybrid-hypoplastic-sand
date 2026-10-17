function [ev,es,e] = inv_e(eps)

% INV_E: Invariants of the strain tensor eps
%        stored as a (6x1) vector in Voigt notation
%        note: shear components doubled
%
%        ev = -tr(eps)   (compresion NEGATIVA; ev > 0 en contraccion; 2026-09-19, rep. 152)
%        q=sqrt((2/3)*tr(e*e))
%


% INV_E: Invariants of the strain tensor epsil stored as a (6x1) vector
%        Returns eps_v=tr(eps), eps_s=sqrt(2/3)|e|,
%        J2e=e_ij*e_ji and J3e=e_ij*e_jk*e_ki
%        model: LaGioia/Nova (1995)
%      
%
fact=1/3;
m=[1;  1;  1;  0;  0;  0];
Minv=[1,  0,  0,   0,   0,   0;
      0,  1,  0,   0,   0,   0;
      0,  0,  1,   0,   0,   0;
      0,  0,  0, 0.5,   0,   0;
      0,  0,  0,   0, 0.5,   0;
      0,  0,  0,   0,   0, 0.5];

ev=-m'*eps;   % 2026-09-19 (rep. 152): ev = -tr(eps), volumetrica ingenieril (contraccion > 0) en compresion NEGATIVA

e=eps+(ev*fact)*m;   % desviador: eps - (tr eps/3) I = eps + (ev/3) I  (2026-09-19)
J2e=e'*(Minv*e);
es=sqrt(2*J2e/3);


