function [Lep,Nep,Hep,z_dot, l_dot, Se, dE,p,FF,dZZ,d_strain2,l_dot2] = hyp_stif(y,parms,deps)
% 
% HYP_STIF: computes hypoplastic tensors L and N
%
% CONVENCION DE SIGNOS (2026-09-19, rep. 152): compresion NEGATIVA en sig y deps, la de
% von Wolffersdorff (1996), Niemunis & Herle (1997) y Liao et al. (2022). p = -tr(sig)/3 > 0
% (inv_s), T_hat = T/tr(T) es invariante, L y N son PARES; sig, deps, T, T_star, delta, z,
% n = de/||de||, z_dot y Hep son IMPARES. El laboratorio (drivers updateModel*, campanas,
% datos, figuras) habla en compresion POSITIVA y mapea en la frontera con espejo_estado.m.
%

M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

I=[1,  0,  0,  0,     0,     0
   0,  1,  0,  0,     0,     0
   0,  0,  1,  0,     0,     0 
   0,  0,  0,  0.5,   0,     0
   0,  0,  0,  0,     0.5,   0
   0,  0,  0,  0,     0,   0.5];

kron_delta=[1,  1,  1,  0,  0,  0];

ny = max(size(y));
sig = y(7:12,1);
DD  = y(1:6,1);
qint = y(13:ny,1);
% Strain_total = qint(2:7,1);

%Ftensor
zint = y(20:25,1);
%Semifluized state
lint = y(26,1);
% 117-B3 / 145 (2026-09-12): dominio de l (Barrero et al. 2020 p.378: l no sale de [0,1]).
% AVISO unico por sesion, no error: un subpaso de prueba del RK puede asomarse fuera y ser
% rechazado por el control de error; abortar ahi cambiaria corridas que hoy terminan.
persistent l_dominio_avisado
if (lint < -1e-9 || lint > 1 + 1e-9) && isempty(l_dominio_avisado)
    l_dominio_avisado = true;
    warning('hyp_stif:l_dominio', 'hyp_stif: l = %.6g fuera de [0,1] (Barrero et al. 2020 p.378). Bajar err_tol.', lint);
end
lint2= y(31,1);

% En   = y(27,1);
% Strain_cum = y(28,1);

%Magnitud del cambio del tensor zint
Zdd = y(29,1);

% add cohesion
sig = sig - kron_delta'*parms(17)/tand(parms(4));  % c is now parms(17). 2026-09-19 (rep. 152): en compresion NEGATIVA la cohesion desplaza sigma en -c/tan(phi) I

e=qint(1);
[p,qdev,z]   = inv_s(sig);
[dev,des,de] = inv_e(deps);
[DDv,DDs,~]  = inv_e(DD);
[Dv,Ds,~]    = inv_e(deps);

T=sig;
T_str=T/(-3*p); % normal stress normalized: T_hat = T/tr(T), tr(T) = -3p en compresion NEGATIVA (2026-09-19, rep. 152). T_hat es PAR: I1, I2, I3, a_1, F y el angulo de Lode no cambian
% Definición de tensores 
T_star = sig+p*kron_delta';         %parte desviatórica  [2026-09-19: sig - (tr sig/3) I = sig + p I]
T_str_star = T_str-kron_delta'/3;   %desviatorica corregida

[p_str,qdev_str,z_str] = inv_s(-T_str);   % 2026-09-19: inv_s recibe tension del MODELO (compresion negativa); T_hat tiene traza +1, se le pasa -T_hat (p_str = 1/3 como antes; salidas no usadas)

TT_str=stress_mult(T_str,T_str);
TT=stress_mult(T_str_star,T_str_star);
TTT=stress_mult(TT,T_str_star);

Tstarnorm = sqrt(T_star'*M*T_star); %norma tensor desviatorico
Tstrstarnorm = sqrt(T_str_star'*M*T_str_star); %norma tensión desviatorica corregida

phic        = parms(4)/180*3.14159265358979;
hs       	= parms(5);
nparam     	= parms(6);
ed0 		= parms(7);
ec0	     	= parms(8);
ei0	   	    = parms(9);
alpha	   	= parms(10);
bphi	   	= parms(11);  % rampa de friccion phi_eff (renombrada de 'beta' 2026-08-08: libera el simbolo del exponente de vW)
% hs_reverse removed from parms vector (not used in model)

tanpsi=sqrt(3.)*Tstrstarnorm;
TTtrace=stress_trace(TT);
if stress_trace(TT)==0
  TTtrace=0.0000000000000001;
end
%% Barotropy, Picknotropy and Volumetric Factors

cos3theta=-sqrt(6.)*stress_trace(TTT)/(TTtrace^1.5);
%Matsuoka Nakai critical state surface criterium
% F1=sqrt((tanpsi*tanpsi/8)+((2-tanpsi*tanpsi)/(2+sqrt(2.)*tanpsi*cos3theta)))-(tanpsi/(2*sqrt(2.)));

% critic friction angle 
% phic = phic;

% phic = phic*(1-bphi) + bphi*phic*y(30,1);
phic = phic*( (1-bphi) + bphi*y(30,1));


%Void ratios definition
expe=-1*(3*p/hs)^nparam;
ei=ei0*exp(expe);
ec=ec0*exp(expe);
ed=ed0*exp(expe);

% smalldiff=0.0001;
% re=(e-ed)/(ec-ed);
% fd=re^alpha;

% if e<=(1+smalldiff)*ed
%   e = ed;
%   disp('e lower than ed');
% end
% 
% if e>=(1+smalldiff)*ei
%   e = ei;
%   disp('e lower than ed');
% end

r = sqrt(T_star'*M*T_star)/p;
% rd   = (ei-parms(32))/(ei-ed);

% I_pe = (ei/e)*exp((ec-e))/((1+r));
% I_pe = 1.05*rd^(rd*alpha)*(ei/e)*exp((ec-e)^3)/((1+r));
% I_pe = rd^(rd*(1-alpha))*(ei/e)^(1+rd)*exp((ec-e))/((1+r)^beta);
% I_pe = 1;
% new expression
I_pe = ((ei/e)*exp((ec-e))/((1+r)))^exp(ec-e); 


%stress invariants Yao et al., 2020
% if qdev > 1e-6
I1 = T_str(1) + T_str(2) + T_str(3);
% CORRECCION E3b (2026-07-02): I2 con convencion estandar (tesis Ap.A, ec.
% Qian et al. 2023 eqs.18-22): I2 = T1*T2+T2*T3+T3*T1 = 0.5*(I1^2 - tr(T^2)).
% La version previa 0.5*(tr(T^2)-I1^2) tenia el SIGNO INVERTIDO y corrompia
% a1 para Omega>0 (phi_TC derivaba, extension no alcanzaba MN). Con Omega=0
% este cambio es matematicamente inerte (a01=1 identicamente).
I2 = 0.5*(I1^2 - T_str'*M*T_str);
% CORRECCION E4-P11 (2026-07-06): I3 = det(T_str) COMPLETO (M-01 auditoria
% externa = P11 reporte 19). La version previa T1*T2*T3 solo es valida SIN
% corte (tensor diagonal); con corte omitia los terminos -T1*T23^2
% -T2*T31^2 -T3*T12^2 +2*T12*T23*T31, corrompiendo a1/a01 para Omega>0
% fuera de ejes principales (I1 e I2 ya eran completos). Orden Voigt
% verificado [11,22,33,12,23,31] (inv_s: t1). INERTE con Omega=0 (a01=1) e
% inerte en estados sin corte (TCU triaxial, plano pi en ejes principales):
% compare esperado sin cambios. Test: audit/E4_P11_test.m
I3 = T_str(1)*(T_str(2)*T_str(3) - T_str(5)^2) ...
   - T_str(4)*(T_str(4)*T_str(3) - T_str(5)*T_str(6)) ...
   + T_str(6)*(T_str(4)*T_str(5) - T_str(2)*T_str(6));
% I3 = T_str(1)*T_str(2)*T_str(3);  % (pre-E4-P11: invalido con corte)

% GUARDAS E3b.2 (2026-07): en el eje hidrostatico (q=0) la formula de a1 es
% 0/0 con limite analitico a1 -> 1 (expansion: num/den -> 4/9, a12 -> 2,
% a13 -> 0). Sin guardas, a0=NaN y a01=1/exp(NaN*Omega)=NaN INCLUSO con
% Omega=0 (NaN*0=NaN), lo que hacia singular a Lep. a0=0 implementa el
% limite exacto; tambien protege sectores con argumento de raiz negativo.
num_a1 = (I1*I2 - I3)*(I1^2 - 3*I2);
den_a1 = I1*I2 - 9*I3;
a0 = 0;                                   % default = limite hidrostatico (F=1)
a0_fallback = 1;   % E4/M-05: 1 = la guarda dejo a0=0 (se refina abajo)
if abs(den_a1) > 1e-14 && (num_a1/den_a1) > 0
    a12 = 3*sqrt(num_a1/den_a1);
    a13 = sqrt(I1^2 - 3*I2);              % = q del tensor normalizado
    if abs(a12 - a13) > 1e-14
        a1 = 2*I1/(a12 - a13);
        if a1 > 0
            a0 = log(a1);
            a0_fallback = 0;              % formula evaluada normalmente
        end
    end
end
% M-05 (2026-07, auditoria externa): el fallback EN el eje hidrostatico
% (q_norm ~ 0) es el limite analitico correcto (E3b.2) y NO se cuenta;
% FUERA del eje era silencioso -> ahora queda registrado en el log
% (log_internal_vars.a0_fallback) para observabilidad. Solo logging:
% no altera ningun resultado. Relevante para Omega>0 en extension.
if a0_fallback == 1 && (I1^2 - 3*I2) <= 1e-12   % q_norm <= 1e-6: eje
    a0_fallback = 0;                             % esperado, no contar
end

% a01= max(real(1/exp(a0*0)),1.0); %valor móvil (0.0,0.5 and 1.0) 
% a01= real(1/exp(a0*0.0)); %valor móvil (0.0,0.5 and 1.0) 


% Generalized failure surface parameter Omega (Yao et al., 2015; Qian et al., 2023)
% Omega=0 -> Drucker-Prager (circular in pi-plane)
% Omega=1 -> Matsuoka-Nakai (curved triangle in pi-plane)
if numel(parms) >= 30
    Omega_param = parms(30);   % <== Omega REAL viene de AQUI (input_data parms(30)). Con canon 42 SIEMPRE esta rama.
else
    Omega_param = 0;    % fallback DP si parms<30 (dead code con canon 42). Antes 0.2 con comentario "DP" erroneo.
end
a01= real(1/exp(a0*Omega_param));







% CORRECCION E3c (2026-07-02): r = Ki/Gi es CONSTANTE MATERIAL.
% Fuente primaria Qian, Wu et al. 2023 (IJNAMG 47:3044, DOI 10.1002/nag.3614):
%   ec.23: a_g = sqrt(3)(3-sin phi_c)/(2 sqrt(2) sin phi_c) * 1/exp(a0*eta)
%          -> el factor exp(a0*eta) entra SOLO en a_g (linea 181 abajo).
%   ec.32: f_v = (3/2)*r - f_a, con "the assumption that Ki/Gi = r is a
%          constant value is adopted".
%   ec.33: f_s = 2(1+e)/(3*lambda*r), sin factor F.
% La version previa r_parms = parms(18)*a01 escalaba r con F(Omega), lo que
% distorsionaba f_s (prop. 1/r) y f_v (prop. r) para Omega>0: la superficie
% real del codigo (E3 Metodo A) no coincidia con el locus declarado
% (Metodo B), max|A-B| ~ 0.4-0.5. Hallazgo P2 del reporte 14, confirmado.
% INERTE con Omega=0 (a01=1 identicamente): tesis y journal no afectados.
% Consistente con tesis cap.3 (F solo en a_g) y journal sec_03 (v_i cte).
r_parms = parms(18); %Material constant Ki/Gi (r_parm is now parms(18))
% r_parms = parms(18)*a01;  % (pre-E3c, incorrecto para Omega>0)


exp1 = nparam*e*(3*p/hs)^nparam;
% exp1 =(3*p/hs)^nparam;
% exp2 = nparam*e*(3*p/(parms(17)*hs))^(nparam);

% CORRECCION E3d (2026-07-02): exponente de a01 en a_g INVERTIDO.
% Cadena de auditoria (doc matlab/audit/E3d_ag_exponente.md):
%  (1) En esta familia L,N vale a_g*M_crit = 3*sqrt(6)/2 exacto (validado
%      E3 con Omega=0 a precision de maquina).
%  (2) Con ag = base*a01 (a01<1 en extension) => M_A = Mc/a01: superficie
%      inflada hacia afuera => alas + NaN de E3 (prediccion rA = rDP^2/rB
%      verificada con error <3% cerca de vertices, 86/98 puntos <8%).
%  (3) El locus declarado (tesis, Metodo B, F_d del journal) es M = Mc*a01
%      (a01 = a1^-Omega, a1>1 en extension = factor de expansion de la
%      tension transformada, Yao 2015). Para producirlo el nucleo requiere
%      a_g = base*exp(+a0*Omega) = base/a01.
% Con Omega=0 es inerte (exp(0)=1). Nota: la ec.23 tipografiada de Qian
% 2023 (1/exp(a0*eta)) contiene la misma trampa de convencion; el ancla es
% la geometria MN + el locus validado + phi_TE(MN)=phi_c.
ag = sqrt(3)/2*(3-sin(phic))/(sqrt(2)*sin(phic))*exp(a0*Omega_param); %E3d: base/a01
% ag = sqrt(3)/2*(3-sin(phic))/(sqrt(2)*sin(phic))*a01; % (pre-E3d, invertido)
% ag = sqrt(3)/2*(3-sin(phic))/(sqrt(2)*sin(phic)); 
% f_a= 1+1/3*ag^2-sqrt(3)/3; %Herle and Kolymbas, 2004
f_a= 1+1/3*ag^2-sqrt(3)/3*ag; %Herle and Kolymbas, 2004
f_v= 3/2*r_parms - f_a;
fd   = (e/ec)^(alpha);

% %Original parameters
% a=sqrt(3.)*(3-sin(phic))/(2*sqrt(2.)*sin(phic));
% NOTA 2026-08-08: el 'beta' de la linea siguiente es el EXPONENTE DE PICNOTROPIA de vW (f_e), otro objeto; NO es parms(11)=bphi.
% fs=(hs/nparam)*(ei/e)^beta*((1+ei)/ei)*(3*p/hs)^(1-nparam)/(3 + a*a - a*sqrt(3.)*((ei0-ed0)/(ec0-ed0))^alpha);

% fs = 2*(1+e)/(3*r_parms)*(1/exp1 - 1/exp2);
% fs = 2*(1+e)/(3*r_parms)*(1/exp_lambda);
fs = 2*(1+e)/(3*1.0*r_parms)*(1/exp1);
% New definition tensors - Ene 2025
Lep = I_pe*fs*(I*(3*p) - f_v*T*kron_delta + ag^2*(T*T')/(3*p));   % 2026-09-19 (rep. 152): T impar => -f_v*T*kron; L es PAR. Autores: L = -I_pe f_s tr(T) [I + f_v T_hat(x)I + a_g^2 T_hat(x)T_hat]
Nep = -I_pe*fs*fd*(ag*(T + T_star));   % 2026-09-19 (rep. 152): T, T_star impares => N es PAR. Ley de los autores: sigma_dot = L:D + N||D|| (Liao 2022 Ec. 1); M_istr suma +rho^chi N(x)delta_hat

% complete von wolferdorf model
% Lep= fs*(F1*F1*I + a*a*T_str*T_str')/stress_trace(TT_str);
% Nep= fs*fd*F1*a*(T_str+T_str_star)/stress_trace(TT_str);

% Wang, Wu, Zhang and Kim (2020)- Overconsilidated clays
% lambda_oed= 0.06;
% kappa     = 0.025;
% yy = (lambda_oed+kappa)/(lambda_oed-kappa);
% c1 = 2*sqrt(3)*ag*(yy-1)/(9*r_parms);
% c2 = 1/3*(sqrt(3)*ag*yy-ag^2-3);
% c3 = ag^2+3*(1-c1);
 
% fs = 3*sqrt(3)*p/(2*ag)*(1/lambda_oed + 1/kappa);
% Lep = (fs)*(I*c1 + c2*T_str*kron_delta + c3*(T_str*T_str'));

%Flow rule m and limit stress condition Y
% qqa = F^2-ag^2*Tstrstarnorm^2/(F^2*(1+c2) + ag^2*sqrt(T_str'*M*T_str)^2);
% m = -T_str_star - qqa*T_str;

% I1= p*3;
% I2= T(1)*T(2)+ T(2)*T(3) + T(3)*T(1);
% I3= T(1)*T(2)*T(3);
% Y = ((yy-1)*(I1*I2/I3-9)+8*tan(phic)^2)/(8*yy*tan(phic)^2);
% fs = 3*sqrt(3)*p/(2*ag)*(1/lambda_oed - 1/kappa);
% Nep = fs*(*m/sqrt(m'*I*m));

% e   = ec;

% using factors by wu
% Lep = (fs)*(I*(3*p_str) + T_str*kron_delta + ag^2*T_str*T_str'/(3*p_str));
% Nep = (fs)*fd*(ag*(T_str + T_str_star));
% e   = ec0*exp(expe);

% using Von Wolferdorf f_s
% Lep = (fs)*I_pe*(C1*I*(3*p_str) + f_v*C2*T_str*kron_delta + C3*T_str*T_str'/(3*p_str));
% Nep = (fs)*I_pe*fd*(C4*(T_str + T_str_star));

% complete wu bahuer
% Lep= fs/I_pe*(a*a*I + T_str*T_str')/stress_trace(TT_str);
% Nep= fs/I_pe*fd*a*(T_str+T_str_star)/stress_trace(TT_str);

global log_internal_vars log_index
log_internal_vars.fs(log_index) = fs;
log_internal_vars.I_pe(log_index) = I_pe;
log_internal_vars.f_v(log_index) = f_v;
log_internal_vars.a01(log_index) = a01;
log_internal_vars.a0(log_index) = a0;
log_internal_vars.a0_fallback(log_index) = a0_fallback;  % M-05 (E4)
log_internal_vars.ag(log_index) = ag;

%% Post liquefaction parameters
cz          = parms(20);  % was 21
zmax        = parms(21);  % was 22
cl          = parms(22);  % was 23
nl          = parms(23);  % was 24
pth         = parms(24);  % was 25
cr          = parms(25);  % was 26
lambda1     = parms(26);  % was 27
lambda2     = parms(27);  % was 28
% lambda3     = parms(32);

Mc = 6*sin(phic)/(3-sin(phic));
% CORRECCION E4-P10 (2026-07-02): Fd debe MULTIPLICAR por a01 = F(Omega).
% El locus declarado de transformacion/falla es q/p = Mc*fd*F (tesis cap.5
% ec. F_d; journal idem; superficie validada E3d: M = Mc*a01). La version
% previa DIVIDIA por a01 => umbral inflado en extension con Omega>0: la
% fabrica z se activaba tarde o no se activaba (C-01 auditoria externa
% 2026-07 = P10 reporte 15). INERTE con Omega=0 (a01=1 identicamente):
% tesis y journal (Omega=0) no afectados. Test: audit/E4_P10_test.m
% ===== SONDA fdDP (2026-08-02): por donde entra el dano de Omega =============
% H2 midio degradacion MONOTONA con Omega acoplado (t18/19/20, Omega 0/0.5/1:
% ext m/e 0.428 -> 0.140 -> 0.094; asim -0.03 -> +0.25 -> +0.29 contra blanco
% -0.38). Hipotesis: el dano no entra por el locus de falla (ag, arriba, que
% conserva a01) sino por AQUI: a01 < 1 en extension ENCOGE la superficie
% Fd = 0 => la extension cruza a dilatancia antes => se corta el tramo
% semifluido de p ~ 0 => la rama extensiva se ARRESTA y el lazo se encoge.
% La sonda SEPARA los dos caminos: falla Matsuoka-Nakai (ag conserva a01) y
% PT Drucker-Prager (Fd sin a01) -- la separacion resistencia/acumulacion de
% Tafili et al. (2024), C&G 171:106314. Con Omega = 0 (produccion) las dos
% ramas son BIT-EXACTAS (a01 = 1 identicamente, y x*1.0 == x en IEEE): E01 y
% S7 ejecutan la rama nueva y son gate real del cambio. H2 lee este flag del
% FUENTE y agrega _fdDP al TAG para no pisar las corridas om1 existentes
% (regla del report 56 6bis). Volver el flag a false restaura la convencion
% declarada E4-P10 (tesis/journal); con Omega = 0 da exactamente lo mismo.
FD_SIN_A01 = false;  % sonda MEDIDA 2026-08-02: trinquete 1-way a extension y densos sin licuar => refutada. false = E4-P10 (produccion)
if FD_SIN_A01
    Fd = (qdev/p)/(Mc*fd)-1;       % PT sin Lode (era la variante comentada)
else
    Fd = (qdev/p)/(Mc*fd*a01)-1;   % E4-P10: produccion declarada
end
% Fd = (qdev/p)/(Mc*fd/a01)-1;  % (pre-E4-P10, incorrecto para Omega>0)
% Fd = (qdev/p)/(Mc*fd*F1)-1;


% Equations Fabric change effect
if sqrt(deps'*I*deps)  ~=  0
        dev_strain_rate_norm = sqrt(de'*I*de);
% CORRECCION E4-P13 (2026-07-04): guarda para 0/0 latente (C-02 auditoria
% externa = P13 reporte 16). Con tasa PURAMENTE volumetrica (deps~=0 pero
% de=0) la normalizacion daba dev_strain_rate = 0/0 = NaN(6,1); el NaN
% alcanzaba z_dot (=> F(20:25)) en estado cizallado con Fd>0, y Ad/Bnew mas
% abajo. Limite adoptado: sin componente desviadora no hay direccion de
% fabrica => dev_strain_rate = 0, LA MISMA convencion que la rama deps==0
% de este mismo if (dev_strain_rate = zeros(6,1)). Nota: el termino
% objetivo zmax*dev_strain_rate es direccional y no tiene limite unico en
% de->0 (misma familia de discontinuidad que P12); la guarda solo hace
% explicita la convencion ya declarada por el codigo. Test: audit/E4_P13_test.m
        if dev_strain_rate_norm > 0
            dev_strain_rate = de/dev_strain_rate_norm;
        else
            dev_strain_rate = zeros(6,1);
        end
        % dev_strain_rate = de/dev_strain_rate_norm;  % (pre-E4-P13: 0/0 => NaN)
        norm_deps = sqrt(deps'*I*deps);        
    if Fd > 0 %dilatation zone implies that z envolves only after phase transformation
%         z_dot = 0*Fd*norm_deps*(zmax*dev_strain_rate - zint);
        z_dot = cz*Fd*norm_deps*(zmax*dev_strain_rate - zint);
    else
        z_dot = zeros(6,1);
    end
        
else
    dev_strain_rate_norm = 0;
    dev_strain_rate = zeros(6,1);
    z_dot = zeros(6,1);
end


%% new state variable
if zint'*deps > 0 && sqrt(zint'*zint) >= Zdd
%     dZZ = 0*sqrt(z_dot'*z_dot);
    dZZ = sqrt(z_dot'*z_dot);
%     dZZ = sqrt(zint'*zint);
else
    dZZ = 0;
end

% FF = (1 - dZZ/zmax);
FF = (1 - Zdd/zmax);
%%
% ===== SONDA g_fl (2026-08-02): Lode SOLO en A_d -- Tafili et al. (2024) ====
% Las tres vias por el TENSOR BASE quedaron medidas y refutadas en H2:
%   Omega acoplado (ag+Fd): arresto, sesgo a compresion (ext 0.43->0.14->0.09
%                           en dosis 0/0.5/1, monotono);
%   fdDP (solo ag, PT DP):  trinquete 1-way a extension en sueltos (asim -1.9
%                           a -3.2, lazos finos DA~0.09) y densos SIN licuar
%                           (t18 121 / t19 37 / t20 25 ramas).
% => la asimetria no se puede pagar por falla/rigidez/PT: un dial, tres precios.
% Tafili et al. (2024), C&G 171:106314, ec.(14): el Lode entra SOLO en el
% amplificador de fabrica,
%     A_d = 1 + <-z:n>/g_fl(theta)
%     g_fl = 2 c_fl/((1+c_fl)-(1-c_fl) cos3theta)   [forma DM04]
%     c_fl = c_mn (3-sin phic)/(3+sin phic),  c_mn = 0.925 (su Fig. 3)
% Compresion: g_fl = 1 (nada cambia). Extension: g_fl = c_fl < 1 y el termino
% de fabrica se amplifica 1/c_fl (~1.55 con 0.925) => mas contraccion => p cae
% mas => tramo semifluido mas largo => MAS acumulacion extensiva. Locus, PT,
% rigidez y monotono INTACTOS ('without spoiling the monotonic behavior').
% IMPLEMENTACION: en trayectoria triaxial el estado es axisimetrico y
% cos3theta = +-1, o sea solo importan los DOS extremos; se interpola con la
% funcion de Lode YA VALIDADA del nucleo (exp(-a0) = M/Mc, cadena
% E3b/E4-P11/E3d) en vez de recomputar cos3theta con otra convencion de
% signos. Extremos EXACTOS de Tafili; forma intermedia MN (irrelevante en
% triaxial; declararlo si algun dia va al UMAT 3D).
% Apagado (0) es BIT-EXACTO: se ejecuta la linea original de produccion.
% Rango util (0, 1.436): mas chico = mas asimetria; 1.436 => sin efecto.
GFL_CMN = 0;   % CERRADO 2026-08-03: dial SIN ventana util (0.925 corto; 0.8 y 0.6 trinquete compresivo, piso 31-66%). 0 = produccion bit-exacta
if zint'*I*dev_strain_rate >= 0
    Ad = 1;
elseif GFL_CMN > 0
    gMN  = real(exp(-a0));                      % M(theta)/Mc del nucleo (1 en TXC)
    gE   = (3 - sin(phic))/(3 + sin(phic));     % = Me/Mc (0.696 con phic = 32.5)
    c_fl = GFL_CMN*gE;                          % Tafili: c_fl = c_mn*(3-sin)/(3+sin)
    g_fl = 1 - (1 - gMN)/(1 - gE)*(1 - c_fl);   % 1 en TXC, c_fl en TXE
    g_fl = min(max(g_fl, min(c_fl, 1)), 1);     % clamp numerico a [c_fl, 1]
    Ad   = 1 + (-zint'*I*dev_strain_rate)/g_fl; % Tafili ec.(14)
else
    Ad = 1 + -zint'*I*dev_strain_rate;          % produccion (Liao 2022)
% Ad = 1;
end


B  = Lep\Nep;
Bp = 1/3*kron_delta'*stress_trace(B); % I*tr(B)/3
Bs = B-Bp;
% ===== R3 (2026-07-16): taper del flujo volumetrico hacia PT — parms(37,38) =====
% h_pt = h0 + (1-h0)*<-Fd>^c37 multiplica SOLO la parte volumetrica Bp (el
% mismo slot quirurgico de A_d, Liao 2022). parms(37)=c37=0 o ausente =>
% h_pt=1 BIT-EXACTO legado. Motivacion (E15/reporte 30): el runaway de
% contractividad a CSR alto ocurre porque tr(B) no decae hasta ~critico
% (locus e-independiente, doc 24 T6 / Fuentes Ap.D); el taper anula el flujo
% volumetrico en el locus e-dependiente Fd=0 (linaje DM04: d ~ A_d(M^d-eta)).
% El piso h0=parms(38) conserva contraccion de reversion y mariposa post-PT
% (taper puro h->0 estanca ru~0.6, medido en port E15). Set unico validado
% en port: c37=1.0, h0=0.25 (con parms(39)=0.008, chimax=13, zmax=9.11).
if numel(parms) >= 37, c37_r3 = parms(37); else, c37_r3 = 0; end
if numel(parms) >= 38, h0_r3  = parms(38); else, h0_r3  = 0; end
% ===== SEMIFLUIDO tent (2026-07-21, reporte 41): parms(42) = tent_mode ==========
% Objetivo: reducir la SOBREESTIMACION de deformacion axial post-licuacion
% (t10: nominal 5x, tent 3x el experimento; figura E25). El taper nominal usa
% <-Fd> (cierra a h0 en la rama dilatante Fd>0 = re-inflacion), suprimiendo la
% dilatancia justo cuando el suelo re-infla p desde ~0 => escape caro. tent usa
% |Fd|: LIBERA la dilatancia en la re-inflacion (h_pt->1 con |Fd| grande) => p
% se recupera con menos corte => excursion menor. NO confundir con liberar en la
% REVERSION (htaper_release, REFUTADO en E24: empeora a 56%).
%   parms(42)=0 o ausente => tent OFF => taper nominal <-Fd> BIT-EXACTO.
%   parms(42)=1 => tent  |Fd|.   parms(42)=2 => tent condicionado a p<pth
%                                (solo regimen semifluido; protege matriz CSR alto).
if numel(parms) >= 42, tent_r41 = parms(42); else, tent_r41 = 0; end
if c37_r3 ~= 0
    if tent_r41 >= 1 && (tent_r41 == 1 || p < pth)
        h_pt = h0_r3 + (1 - h0_r3)*min(abs(Fd), 1)^c37_r3;   % tent: libera re-inflacion
    else
        h_pt = h0_r3 + (1 - h0_r3)*max(-Fd, 0)^c37_r3;       % nominal <-Fd> (bit-exacto)
    end
else
    h_pt = 1;
end
Bnew = Ad*h_pt*Bp + Bs;
% Bnew = Ad*Bp + Bs;   % (legado pre-R3; identico con parms(37)=0)

% CORRECCION E4-P12 (2026-07-06): gate isocorico con tolerancia RELATIVA.
% Antes 'if Dv == 0' (igualdad exacta en float): la correccion de Liao
% (escala Ad solo en la parte volumetrica de B) se activaba SOLO con
% tr(deps) = 0 bit-exacto. En control mixto deps sale del NR (tol 1e-10) y
% tr(deps) queda como residuo != 0 => el gate se apagaba de forma
% DISCONTINUA ante perturbaciones infinitesimales (M-02 auditoria externa
% = P12). Criterio: isocoria numerica |Dv| <= tol*||deps|| (adimensional).
% tol = 1e-8 separa residuos numericos (ratios <= 8.5e-11 medidos en
% TCU_baseline_E4.csv) de trayectorias volumetricas reales (ratio O(1)).
% Caso deps = 0: 0 <= 0 activa el gate, identico a la version previa.
% Verificacion: E01 + driver 'compare' (esperado max|dif| <= 1e-10; en TCU
% Dv = 0 exacto => gate ya estaba activo => sin cambios).
tol_Dv_gate = 1e-8;
if abs(Dv) <= tol_Dv_gate*sqrt(deps'*I*deps)
Nep  =  Lep*Bnew;  %signo negativo en -Lep (Ver artículo de Liao et al., 2022)
end
% if Dv == 0   % (pre-E4-P12: igualdad exacta, gate discontinuo)

% Dc = -B;

%Semifluized state
pr = p/(pth);
fl = 0.01;                                   % f_l default con cita (Liao 2022; rep. 141 sec. 2)
if numel(parms) >= 43 && parms(43) > 0       % 145 / H4 (2026-09-12): f_l barrible por parms(43); ausente o 0 => 0.01 BIT-EXACTO
    fl = parms(43);
end
% fl = 0;

% ===== E4-P14 DECIDIDO (2026-07, opcion a): comportamiento INTENCIONAL =====
% En pr = 1 hay salto de Se: rama pr<=1 da Se -> 1 + fl = 1.01; rama pr>1
% da Se = 1. Se DECIDIO conservar el codigo tal cual (validado, calibrado,
% resultados publicados) y corregir en E5 la DECLARACION de tesis/journal
% (deben decir: Se = pr^(lambda*lint) + fl para pr <= 1; Se = 1 para
% pr > 1). La "forma acotada" de la 2a auditoria externa RECALIBRA y fue
% descartada. Nota: pr >= 0 garantizado por clamp p >= 1e-4 (inv_s, P17).
if pr <= 1
    % A22 2026-08-07: termino (p_in/p_inr)^a de Barrero Ec.(5), que Liao omitio.
    % Barrero et al. 2020, IJNAMG 44:371-388. El factor va DENTRO del corchete
    % de <L>: multiplica la generacion, NO el retorno a cero (cr). Verificado.
    % parms(31)=a (0 => APAGADO y BIT-EXACTO, rama else) | parms(32)=p_inr kPa
    % | parms(33)=p_in, que escribe updateModel_stress_ES en cada reversion.
    if parms(31) ~= 0
        if numel(parms) < 33 || parms(32) <= 0 || parms(33) <= 0
            error('hyp_stif:A22_sin_pin', ...
                'A22: a=%g activo pero p_inr/p_in sin inicializar (falta la Edicion 2).', parms(31));
        end
        fpin_a22 = (parms(33)/parms(32))^parms(31);
        l_gen = fpin_a22*abs(des)*( cl*(1-pr )*(1-lint)^nl );   % 117-B2: solo la generacion; el retorno va fuera de la compuerta
    else
    l_gen   = abs(des)*( cl*(1-pr )*(1-lint)^nl );   % 117-B2: solo la generacion
    end
% l_dot   = 0.0*abs(des)*( 0*(1-pr )*(1-lint)^0 ) - 0*lint*abs(dev);
l_dot2  = 0*abs(Nep'*deps)*( 0.3*cl*(1-pr )*(1-lint2)^nl ) - cr*lint*abs(dev);
    
    
lambda= lambda1*exp(lambda2*((ei-e)/(ei-ed)));

%     lambda = lambda1*exp(((ei-e)/(ei-ed)));  %Corregir se fuera necesario
    % bbb = 0.9;
    bbb = 1;

    Se    = (1 - (1-pr))^(lambda*bbb*lint) + fl;
%     Se    = 1;

else
    
    l_gen = 0;   % 117-B2
    l_dot2= 0;
    Se    = 1;

end
% 117-B2 / 145 (2026-09-12): retorno a cero de l FUERA de <1-pr> (Barrero et al. 2020 ec. 5 y p.377;
% Liao et al. 2022 ec. 10). Bit-exacto en la campana: no drenado => tr(deps) = 0 exacto; drenado y
% edometrico => p > p_th siempre => l = 0 => el termino vale 0.
l_dot = l_gen - cr*lint*abs(dev);

% if DDs > eq_lim
% Br= En/(En+eq_lim);
if Dv > 0
    Dv = Dv;
else
    Dv = 0;
end

dE= Dv*p/(3*p)+qdev*Ds/(3*p); %particle breakage
d_strain2 = dev_strain_rate_norm;
% d_strain2 = parms(20) *(1 - Strain_total2)*dev_strain_rate_norm;
% e = e*exp(Br*e_br);

% end
% dE = 0;
Hep = zeros(6,1);
Hep(1:3) = (1+e);   % 2026-09-19 (rep. 152): e_dot = (1+e) tr(D) en compresion NEGATIVA (contraccion: tr D < 0 => e baja); antes -(1+e) con tr D > 0




% log_internal_vars.z_dot(log_index) = z_dot;
log_internal_vars.eps_a(log_index) = y(3); 
% log_internal_vars.Se(log_index) = Se;
log_internal_vars.FF(log_index) = FF;
log_internal_vars.dZZ(log_index) = dZZ;
log_internal_vars.qdev(log_index) = qdev;
log_internal_vars.p(log_index) = p;
log_internal_vars.pr(log_index) = pr;
% log_internal_vars.l_dot(log_index) = l_dot;
% log_internal_vars.lint(log_index) = lint;
%also
log_internal_vars.Fd(log_index) = Fd;
log_internal_vars.dev_strain_rate_norm(log_index) = dev_strain_rate_norm;
% ===== F1 (2026-07-26): observabilidad de la FABRICA — SOLO LOGGING =====
% Ad y ||z|| no se registraban, asi que la unica forma de saber cuanto trabaja el
% mecanismo era ablacionarlo y volver a correr — caro, porque sin fabrica el
% modelo licua mucho mas tarde y el early-stop no llega a dispararse. Con estas
% dos lineas basta UNA corrida de referencia para medirlo: Ad = 1+<-z:n> es el
% factor por el que se amplifica el flujo VOLUMETRICO (Bnew = Ad*h_pt*Bp + Bs),
% y ||z||/z_max dice cuanto se desarrollo la fabrica.
% NO altera ningun calculo: son asignaciones a la struct global de logging.
log_internal_vars.Ad(log_index) = Ad;
log_internal_vars.z_norm(log_index) = sqrt(zint'*zint);
log_index = log_index + 1;
