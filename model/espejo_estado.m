function y = espejo_estado(y)
% ESPEJO_ESTADO  Mapa exacto entre la convencion del MODELO (compresion NEGATIVA: von
%   Wolffersdorff 1996, Niemunis & Herle 1997, Liao et al. 2022) y la del LABORATORIO
%   (compresion POSITIVA: datos, drivers updateModel*, campanas y figuras). Es una
%   involucion: se aplica igual en los dos sentidos.
%
%   Estado y (31 x 1, o 31 x N por columnas):
%     impares (cambian de signo): 1:6 eps, 7:12 sig, 14:19 delta (deformacion
%                                 intergranular), 20:25 z (fabrica)
%     pares   (no cambian):       13 e, 26 l, 27 dE, 28 eps_cyc, 29 Zdd, 30 eps_mon, 31
%   Se escribe 0 - y y no -y para que las ranuras nulas queden en +0 (nunca -0); para
%   todo numero distinto de cero las dos formas son la misma negacion exacta (IEEE 754).
%
%   2026-09-19, rep. 152 (04_Auditoria/Audit_report/152_MAPA_PARIDAD_CONVENCION_2026-09-19.md)
IMPAR = [1:12, 14:25];
y(IMPAR, :) = 0 - y(IMPAR, :);
end
