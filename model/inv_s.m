function [p,q,z] = inv_s(sigma)

% INV_S: Invariants of the stress tensor sigma
%        stored as a (6x1) vector in Voigt notation 
%
%        p = -tr(sig)/3   (compresion NEGATIVA: p > 0; 2026-09-19, rep. 152)
%        q=sqrt((3/2)*tr(s*s))
%        z = sin(3*theta), theta = Lode angle
%

persistent p_clamp_avisado   % E4-P17: aviso unico del clamp de p
fact=-1/3;   % 2026-09-19 (rep. 152): compresion NEGATIVA, p = -tr(sigma)/3 > 0. El mismo fact en temp1 deja z = sin(3 theta) como antes (+1 en compresion triaxial)
m=[1;  1;  1;  0;  0;  0];
M=[1,  0,  0,  0,  0,  0
   0,  1,  0,  0,  0,  0
   0,  0,  1,  0,  0,  0
   0,  0,  0,  2,  0,  0
   0,  0,  0,  0,  2,  0
   0,  0,  0,  0,  0,  2];

p=fact*m'*sigma;
% ===== E4-P17 (2026-07): clamp de presion media DOCUMENTADO =====
% Politica: p se acota inferiormente a p_tol = 1e-4 kPa. Consecuencias:
% (a) no hay division por cero aguas abajo (T/(3p), q/p, p/pth);
% (b) pr = p/pth nunca es negativo => la potencia Se = pr^lambda es real;
% (c) cerca de licuacion (p -> 0) el estado se evalua en p = p_tol
%     (mal condicionamiento posible de Lep\Nep: vigilar, no error).
% Antes era SILENCIOSO (hallazgo C-03/P17): ahora avisa UNA vez por
% sesion la primera vez que actua. El clamp NO cambia: resultados
% identicos. Documentar tambien en tesis Ap.B (E5).
if p < 0.0001
    if isempty(p_clamp_avisado)
        warning('inv_s:p_clamp', ...
            'inv_s: p = %.3e kPa < 1e-4 -> clamp a 1e-4 kPa (aviso unico por sesion; E4-P17).', p);
        p_clamp_avisado = true;
    end
    p = 0.0001;
end
s=sigma+p*m;   % desviador: sigma - (tr sigma/3) I = sigma + p I  (2026-09-19)
sn=s'*(M*s);

q=sqrt(3*sn/2);

% components of s*s as a (6x1) vector t1

t1=[s(1)*s(1)+s(4)*s(4)+s(6)*s(6);
    s(4)*s(4)+s(2)*s(2)+s(5)*s(5);
    s(6)*s(6)+s(5)*s(5)+s(3)*s(3);
    2*(s(1)*s(4)+s(4)*s(2)+s(6)*s(5));
    2*(s(4)*s(6)+s(2)*s(5)+s(5)*s(3));
    2*(s(6)*s(1)+s(5)*s(4)+s(3)*s(6))];

% Lode angle

if q==0.0
    z = 0;
else
    temp1=fact*t1'*s;
    z=27.0*temp1/(2.0*q*q*q);
    if abs(z) > 1.0
        z=z/abs(z);
    end
end

