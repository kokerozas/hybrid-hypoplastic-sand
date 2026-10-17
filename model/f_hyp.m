function F = f_hyp(y,S,E,V,parms)

% F_EP: computes RHS of the evolution ODE for plastic processes
%

M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

M2=[1,  0,  0,  0,    0,    0
   0,  1,  0,  0,    0,    0
   0,  0,  1,  0,    0,    0 
   0,  0,  0,  0.5,  0,    0
   0,  0,  0,  0,    0.5,  0
   0,  0,  0,  0,    0,    0.5];
kron_delta = [1;1;1;0;0;0];

mR	   	= parms(12);
ny 		= max(size(y));
qint 		= y(13:ny,1);
istrain		= qint(2:7,1);

ny = max(size(y));
F = zeros(ny,1);

[Lep] = hyp_stif(y,parms,[0,0,0,0,0,0]');
Hdel = zeros(6,6);

Maux=Lep*mR;
A    = S*Maux+E;
depsaux = A\V;

[Lep,Nep,Hep,Ftensor_dot,H, Se, dE,p,FF,dZZ,~,ldot2] = hyp_stif(y,parms,depsaux);

% ===== R3b (2026-07-16): chimax por ESTADO — parms(40)=k_chi, parms(41)=psi_ref =====
% chimax_eff = parms(29)*exp(k_chi*(psi - psi_ref)), psi = e - e_c(p) local.
% Resuelve el residuo de densidad-suelta (reporte 31): el techo del ramp
% chi(eps_cyc) es mayor para estados sueltos (mas proteccion en ensayos
% largos), constante para densos. UNA constante nueva (k_chi); psi_ref es el
% estado de anclaje DECLARADO (serie densa TCUI18-20, -0.181), no libre.
% parms(40)=0 o ausente => parms_mi=parms BIT-EXACTO. Sin cambio de firma:
% se pasa una copia local de parms con (29) reescalado a M_istr.
% Port (kchi=3.6): t8 0.76 | t9 0.96 | t10 0.89 | t18 1.03 | t19 1.00 | t20 0.99.
parms_mi = parms;
if numel(parms) >= 40 && parms(40) ~= 0
    psi_r3b = y(13,1) - parms(8)*exp(-(3*p/parms(5))^parms(6));
    if numel(parms) >= 41, psiref_r3b = parms(41); else, psiref_r3b = -0.181; end
    parms_mi(29) = parms(29)*exp(parms(40)*(psi_r3b - psiref_r3b));
end

% B = -Lep\Nep;
% norm_deps = sqrt(depsaux'*M2*depsaux);
% norm_B    = sqrt(B'*M*B);
% norm_Bdeps= sqrt((B'*M2*depsaux)^2);
% cos_w= norm_Bdeps/(norm_deps*norm_B);
% % 
% Nep = Nep*cos_w^2;

% % --- Toggle plástico global (ej.: parms(34)=0 apaga)
% plastic_on = (numel(parms) >= 34) && (parms(34) ~= 0);
% if ~plastic_on
%     Nep(:,:) = 0;   % corta la rama “plástica” hipoplástica
% end

%% compute A matrix and solve for deps

Maux=Lep*mR;
A    = S*Maux+E;
depsaux = A\V;

[Mep,Hdel,dStrain,dStrain2]=M_istr(Lep,Nep,depsaux,istrain,parms_mi,y(28,1),FF,y(30,1));
Mep = Mep*(1/3*(1-Se)*(kron_delta.*kron_delta') + Se*eye(6));   % 145 (2026-09-12): la identidad de 4o orden es eye(6) en Voigt ingenieril; con M2 el bloque de corte recibia Se/2


A    = S*Mep+E;
depsMNR = A\V;
depsinit = depsMNR;
sigpresc = Mep*depsMNR;

[Mep,Hdel,dStrain,dStrain2]=M_istr(Lep,Nep,depsMNR,istrain,parms_mi,y(28,1),FF,y(30,1));
Mep = Mep*(1/3*(1-Se)*(kron_delta.*kron_delta') + Se*eye(6));   % 145 (2026-09-12): la identidad de 4o orden es eye(6) en Voigt ingenieril; con M2 el bloque de corte recibia Se/2

%% Modified Newton-Raphson interations for correction of deps
stresscontrol=zeros(6,1);
for i = 1:6
  for j = 1:6
    if  (S(i,j) > 0)
      stresscontrol(j)=1;
    end
  end
end

isstresscont=sqrt(stresscontrol'*M*stresscontrol);

if isstresscont > 0
  tol_mnriter  = parms(3);
  itererror=1000;
  max_iter = 50;
  iter = 0;
  
  while (itererror > tol_mnriter) && (iter < max_iter)
        iter = iter + 1;
        sigest= Mep*depsMNR;
        unbalsig=sigpresc-sigest;
 
        for i = 1:6
	  if  (stresscontrol(i) == 0)
	    depsMNR(i)=depsinit(i);
	    unbalsig(i)=0;
	  end
        end

        depsMNR=depsMNR + Mep\unbalsig;
        itererror=sqrt(unbalsig'*M*unbalsig);

        [Mep,Hdel,dStrain,dStrain2]=M_istr(Lep,Nep,depsMNR,istrain,parms_mi,y(28,1),FF,y(30,1));
        Mep = Mep*(1/3*(1-Se)*(kron_delta.*kron_delta') + Se*eye(6));   % 145 (2026-09-12): la identidad de 4o orden es eye(6) en Voigt ingenieril; con M2 el bloque de corte recibia Se/2

  end
  if iter >= max_iter
        warning('Modified Newton-Raphson reached maximum iterations (%d) with error %e', max_iter, itererror);
  end
end
%% F vector

F(1:6,1)   = depsMNR;
F(7:12,1)  = Mep*depsMNR;
de= Hep'*depsMNR;
distrain=Hdel'*depsMNR;
F(13,1) = de;
F(14:19,1) = distrain;
F(20:25,1) = Ftensor_dot;
F(26,1) = H;
F(27,1) = dE;
F(28,1) = dStrain;
F(29,1) = dZZ;
F(30,1) = dStrain2;
% ===== E4 paso 3 (2026-07-06, hallazgo C-04/P16): l_dot2 = 0 =====
% y(31)=lint2 era una integracion MUERTA: ningun modulo consume lint2 y la
% norma de error del RK la excluye (hypo_update: active_idx=[20:26,28,30]).
% Se anula para eliminar la deriva inerte. Verificacion: E01 (sin cambios)
% + driver 'compare' (max|dif| <= 1e-10 vs TCU_baseline_E4.csv).
% Linea original:
% F(31,1) = ldot2;%
F(31,1) = 0;



