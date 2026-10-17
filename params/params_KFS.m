function [parms, y0, TESTS] = params_KFS(set)
%PARAMS_KFS  Karlsruhe fine sand (KFS): material constants, initial internal state and published tests.
%
%   [parms, y0, TESTS] = params_KFS('cyclic')     closed cyclic set: undrained cyclic tests (TCU)
%   [parms, y0, TESTS] = params_KFS('monotonic')  monotonic closure: TMD, TMU and OED
%
%   parms  42 slots (column). Slot 43 (f_l) is optional and absent: the core uses f_l = 0.01.
%   y0     31 state variables in the LABORATORY convention (compression positive). The drivers in
%          model/ map it to the model convention (compression negative) with model/espejo_estado.m.
%   TESTS  published undrained cyclic tests (Wichtmann & Triantafyllidis 2016) with the initial
%          state, amplitude and cycle budget used in the thesis.
%
% The two sets differ only in slots 2 (err_tol), 22 (c_l), 23 (n_l), 26 (lambda_1) and 27 (lambda_2).
% Slots 22-27 belong to the semifluidized state and stay inactive while p' > p_th, which is the case in
% the monotonic tests.
%
% parms(34) is ETA, the rate of the cyclic memory y(28), and parms(35) is KAPPA, the rate of the
% monotonic memory y(30). Earlier versions named them kappa_rate and eta_rate (inverted); model/M_istr.m
% reads them as eta_rate_cyc and kappa_rate_mon. The full map of slots is in docs/NOMENCLATURE.md.

if nargin < 1 || isempty(set), set = 'cyclic'; end
assert(any(strcmp(set, {'cyclic', 'monotonic'})), 'params_KFS: set must be ''cyclic'' or ''monotonic''.');

%% integration
max_ksub = 1000;      % maximum number of substeps per step
err_tol  = 1.0e-3;    % normalized error tolerance of the adaptive Runge-Kutta 3(2)
tol_NR   = 1.0e-10;   % tolerance of the Newton iteration under mixed or stress control

%% hypoplasticity (von Wolffersdorff 1996) with the friction ramp
phi_c = 33.11;    h_s = 4e6;    n = 0.27;                 % critical friction angle [deg], granular hardness [kPa], exponent
e_d0  = 0.677;    e_c0 = 1.054;  e_i0 = 1.212;            % limit void ratios at p = 0
alpha = 0.8311;   b_phi = 0.1114;                          % pycnotropy exponent; ramp of the mobilized friction with eps_mon
r_KG  = 1.0902;   c = 0;                                   % ratio K_i/G_i; cohesion (off)

%% intergranular strain (Niemunis & Herle 1997) and cyclic memory
m_R = 2.2;   m_T = 1.1;   R = 1.0e-4;   beta_r = 0.15;   chi_0 = 5.5;   chi_max = 10;

%% fabric (Liao et al. 2022) and semifluidized state (Barrero et al. 2020, form of Liao et al. 2022)
c_z = 500;   z_max = 9.1084;
c_l = 24;    n_l = 16;    p_th = 10;    c_r = 1000;    lambda_1 = 0.627986;    lambda_2 = 1.5;

%% memories: eta (rate of eps_cyc) and kappa (rate of eps_mon)
eta = 20;    kappa = 9.9575;

%% switches and retired slots (kept at the values of the thesis)
Omega = 0;                      % failure surface: 0 = the one of the thesis
a = 0;   p_inr = 0;   p_in = 0; % pressure factor of Barrero et al. (2020): must stay 0 in this release
y_h = 1;   c_h = 0;   h_0 = 0.25;   g = 0;   k_chi = 0;   psi_ref = -0.218;   tent = 2;
q_br = 9.9575;   e0_doc = 0.798; % not read by the core (documentary)

if strcmp(set, 'monotonic')     % monotonic closure (before the cyclic closure)
    err_tol  = 1.0e-4;
    c_l = 7.4;   n_l = 8;   lambda_1 = 0.06369;   lambda_2 = 4.5;
end

assert(y_h == 1 && g == 0 && k_chi == 0, ['params_KFS: y_h = 1, g = 0 and k_chi = 0 give the cyclic ' ...
       'memory of the thesis and the manuscript (detector y_delta, docs/NOMENCLATURE.md).']);
assert(a == 0, 'params_KFS: a must be 0 (p_in at reversals is not implemented in model/updateModel_stress.m).');
parms = [max_ksub; err_tol; tol_NR; phi_c; h_s; n; e_d0; e_c0; e_i0; alpha; b_phi; ...   %  1-11
         m_R; m_T; R; beta_r; chi_0; c; r_KG; q_br; ...                                    % 12-19
         c_z; z_max; c_l; n_l; p_th; c_r; lambda_1; lambda_2; e0_doc; chi_max; Omega; ...  % 20-30
         a; p_inr; p_in; eta; kappa; y_h; c_h; h_0; g; k_chi; psi_ref; tent];              % 31-42

%% initial state (laboratory convention); stress and void ratio are set by init_state.m
eps0   = zeros(6, 1);
sig0   = [100; 100; 100; 0; 0; 0];
e0     = 0.798;
delta0 = -0.82*R/sqrt(3)*[1; 1; 1; 0; 0; 0];              % intergranular strain, isotropic, ||delta(0)|| = 0.82 R
z0     = 0.5*[-1/sqrt(6); -1/sqrt(6); 2/sqrt(6); 0; 0; 0];  % fabric, ||z(0)|| = 0.5 along triaxial compression
y0     = [eps0; sig0; e0; delta0; z0; zeros(6, 1)];        % 26 l, 27 -, 28 eps_cyc, 29 -, 30 eps_mon, 31 -

%% published undrained cyclic tests used in the thesis
%   q_amp = CSR*p0; ncyc = cycles integrated in the thesis run; N_test = cycles to initial
%   liquefaction reported by Wichtmann & Triantafyllidis (2016)
%            id   CSR   e0     Dr   p0   ncyc  N_test
T = [         8   0.20  0.821  62   100   304   249
              9   0.25  0.798  68   100    80   100
             10   0.30  0.825  61   100    17    15
             18   0.30  0.759  78   100    60    54
             19   0.40  0.761  78   100    26    15
             20   0.50  0.755  79   100    12     6 ];
TESTS = struct('id', T(:,1), 'CSR', T(:,2), 'e0', T(:,3), 'Dr', T(:,4), 'p0', T(:,5), ...
               'q_amp', T(:,2).*T(:,5), 'ncyc', T(:,6), 'N_test', T(:,7));
end
