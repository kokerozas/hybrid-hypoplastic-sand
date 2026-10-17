 function [Mep,Hdel,d_strain,d_strain2] = M_istr(Lep,Nep,deps,istrain,parms,Strain_total,FF,Strain_total2)

% computes matrix M of the intergranular strain concept by Niemunis and Herle (1997)
% CONVENCION (2026-09-19, rep. 152): compresion NEGATIVA. Lep, Nep y Mep son PARES; deps y
% istrain (delta) IMPARES. La mezcla lleva +rho^chi N(x)delta_hat como en la Ec. 8 de N&H 1997.

M2=[1,  0,  0,  0,    0,    0
   0,  1,  0,  0,    0,    0
   0,  0,  1,  0,    0,    0 
   0,  0,  0,  0.5,  0,    0
   0,  0,  0,  0,    0.5,  0
   0,  0,  0,  0,    0,    0.5];

I=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  1,  0,  0
   0,  0,  0,  0,  1,  0
   0,  0,  0,  0,  0,  1];

mR= parms(12);
mT= parms(13);
Rparam= parms(14) ;
betar= parms(15);
chimax=parms(29);  % was parms(33) before vector compaction

%Poblete 2016
chi    = (parms(16) + (Strain_total *(chimax - parms(16))));
% chi    = (parms(16) + (Strain_total *(12 - parms(16))));


istrain_norm=sqrt(istrain'*M2*istrain);
TINY=1.e-20;
istrain_dir=zeros(6,1);

if  (istrain_norm>TINY)
     istrain_dir=istrain/istrain_norm;
end

istraindir_D=istrain_dir'*M2*deps;
rho=istrain_norm/Rparam;

if rho > 1
    rho = 1;
end

LijmnSmn=Lep*istrain_dir;
LijmnSmnSkl=LijmnSmn*istrain_dir'*M2;   % 145 (2026-09-12): la contraccion hd:deps lleva M2, la misma metrica que istraindir_D (L38)
NijSkl=Nep*istrain_dir'*M2;    % dyadic product properties; 145: idem, hd:deps con M2

Mep=(mT*rho^chi+mR*(1-rho^chi))*Lep;
% Mep=mT*Lep; %monotonic condition

% ===== R3 (2026-07-16): piso de sub-movilizacion en el slot N — parms(39) ====
% rcN = max(rho^chi, g39) SOLO en el termino -rho^chi*Nep*delta_hat de la
% mezcla: acumulacion pequenia pero no nula dentro de la burbuja IS (deficit
% de shakedown a CSR bajo, E13/E15; fisica de umbral gamma_t). La rigidez L
% (mR/mT/rho^chi) queda INTACTA; monotonico inerte (rho->1 => rcN=rho^chi).
% parms(39)=0 o ausente => BIT-EXACTO legado. Set unico port: g39=0.008.
if numel(parms) >= 39, g39_r3 = parms(39); else, g39_r3 = 0; end
rcN_r3 = max(rho^chi, g39_r3);
if(istraindir_D>0)
  Mep=Mep+(1-mT)*rho^chi*LijmnSmnSkl+rcN_r3*NijSkl;   % 2026-09-19 (rep. 152): +rho^chi N(x)delta_hat (N&H 1997 Ec. 8; Liao 2022 Ec. 2). N PAR, delta_hat IMPAR => NijSkl cambia de signo y el termino queda identico
% Mep=Mep+(1-mT)*rho^chi*LijmnSmnSkl+rho^chi*NijSkl;  % (legado; identico con parms(39)=0)
  Hdel=I-rho^betar*M2*istrain_dir*istrain_dir';
else
  Mep=Mep+(mR-mT)*rho^chi*LijmnSmnSkl;
  Hdel=I;
end



dev_strain_rate_norm = sqrt(deps'*M2*deps);   % 145 / D3 (2026-09-12): norma de Frobenius en Voigt ingenieril; con I = eye(6) los cortantes pesaban doble
% ===== P21 (2026-07): selector de definicion de y_h — adimensionalizacion =====
% parms(36) = yh_mode.  Default 0 (o numel(parms)<36) = LEGADO BIT-EXACTO.
%   0 -> yh = (rho^chi * Nep)' * deps                      [proyeccion sin normalizar; tesis v1, P21]
%   1 -> yh = rho^chi * <delta_hat : deps> / ||deps||_M2   [Variante A, ISA fiel]
% Doc 24 §2: la "N" de Fuentes (2020) Ec.29 es la DIRECCION de flujo
% intergranular delta_hat (= istrain_dir), NO el tensor N hipoplastico (Nep).
% Variante A traduce fielmente Ec.29 al marco Niemunis-Herle (1997):
%   <.> = Macaulay = max(.,0) -> yh=0 exacto en reversion (delta_hat opuesto a deps)
%   / ||deps||_M2             -> coseno: acotada [0,1], paso- y tension-invariante.
% istraindir_D (= delta_hat:deps, metrica M2) ya esta calculado arriba (L38).
% Un unico yh alimenta el log y d_strain: elimina el duplicado silencioso
% (antes la proyeccion se recomputaba inline dentro de d_strain).
if numel(parms) >= 36, yh_mode = parms(36); else, yh_mode = 0; end
norm_deps_M2 = sqrt(deps'*M2*deps);
if yh_mode == 1 && norm_deps_M2 > TINY
    yh = rho^chi * max(istraindir_D, 0) / norm_deps_M2;   % Variante A (ISA fiel)
else
    yh = -(rho^chi*Nep)'*deps;                              % legado exacto (P21); 2026-09-19: N par, deps impar => el signo conserva y_h
end
%deps := incremento de deformación aplicado

% ===== E4-P6 (2026-07): nombres VERACES para las tasas de memoria =====
% Historico: los nombres locales estaban INVERTIDOS respecto de la tesis:
%   parms(34) (ex "kappa_rate") = 20 es eta   (tesis cap.5: tasa de e_cyc;
%                                 alimenta d_strain  -> y(28) "Strain_total")
%   parms(35) (ex "eta_rate")   = 15 es kappa (tesis cap.3: tasa de e_mon;
%                                 alimenta d_strain2 -> y(30) "Strain_total2")
% El MAPEO NUMERICO parms(34)->d_strain, parms(35)->d_strain2 se conserva
% EXACTO (solo se renombran variables locales): resultados identicos.
% Config de tesis: numel(parms)=30 => rigen los defaults 20/15.
if numel(parms) >= 34
    eta_rate_cyc = parms(34);      % eta (tasa e_cyc)   [ex kappa_rate]
else
    eta_rate_cyc = 20;
end

if numel(parms) >= 35
    kappa_rate_mon = parms(35);    % kappa (tasa e_mon) [ex eta_rate]
else
    kappa_rate_mon = 15;
end

% d_strain = 5.0*(1 - yh - Strain_total) * dev_strain_rate_norm;
d_strain = eta_rate_cyc*(1 - yh - Strain_total) *dev_strain_rate_norm;   % P21+F3: eta_rate_cyc=parms(34) (bit-exacto: =20 si numel(parms)<34)
% FASE 0 (2026-07-26): SLOT CONECTADO. Hasta ahora esta linea llevaba un LITERAL
% (15 historicamente, 9.95 desde el 2026-07-26) mientras kappa_rate_mon se leia de
% parms(35) en la linea de arriba y se DESCARTABA: kappa figuraba en la tabla de
% parametros de la tesis pero el codigo no lo leia, y cambiar el slot no hacia nada.
% Ahora rige parms(35). Valor adoptado: 9.9575, que es el calibrado sobre el peak de
% resistencia en triaxiales monotonos drenados (estaba aparcado en el slot 19 bajo el
% nombre "q_br", que el nucleo no usa).
%
% OJO — kappa NO es exclusivo del regimen monotono. Strain_total2 = y(30) modula el
% angulo de friccion critico en hyp_stif L89:  phic = phic*((1-beta) + beta*y(30)),
% con beta = parms(11) = 0.1114 => phi_c va de 28.88 deg (y30=0) a 32.5 deg (y30=1),
% o sea Mc varia un 14%. Y Mc define F_d, que gobierna la fabrica, el taper y tent.
% Como y30 = 1 - exp(-kappa*s), bajar kappa retrasa la saturacion: el efecto es
% APRECIABLE en los ensayos cortos (t19, t20) y NULO en los largos (t8, t9), donde
% y30 satura igual. Verificacion: audit/S4_kappa_check.m
d_strain2 = kappa_rate_mon *(1 - Strain_total2)*dev_strain_rate_norm;

% if Strain_total < 0
%     d_strain = 0;
% end




global log_internal_vars log_index
log_internal_vars.d_strain(log_index) = d_strain;
log_internal_vars.Strain_total(log_index) = Strain_total;
log_internal_vars.Strain_total2(log_index) = Strain_total2;
log_internal_vars.d_strain2(log_index) = d_strain2;
log_internal_vars.chi(log_index) = chi;
log_internal_vars.yh(log_index) = yh;
log_internal_vars.rho(log_index) = rho;                       % R2/E14b: rho CAPADO (el que entra a yh) — canal 1 (movilizacion)
log_internal_vars.rho_raw(log_index) = istrain_norm/Rparam;  % R2/E14b: rho SIN cap (>1 => satura); solo logging, inerte
log_internal_vars.yh_mode(log_index) = yh_mode;   % P21: cada CSV declara que definicion de y_h corrio (inerte)
% disp(['Chi actual: ', num2str(chi), ' | FF: ', num2str(FF), ' | Strain_total: ', num2str(Strain_total), ' | d_strain: ', num2str(d_strain), ' | d_strain2: ', num2str(d_strain2)]);

end
