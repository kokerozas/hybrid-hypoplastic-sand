function [parms, y0, TESTS] = params_BBS(set)
%PARAMS_BBS  Bio-Bio sand (BBS): material constants, initial internal state and published tests.
%
%   [parms, y0, TESTS] = params_BBS('cyclic')    closed set of the 11 m layer (undrained cyclic tests
%                                                of Saldana Sotelo 2023): the set of the thesis for BBS
%   [parms, y0, TESTS] = params_BBS('monotonic') set of the drained tests of another Bio-Bio sample
%                                                (Cuitino et al. 2015; Brown et al. 2024), used in Fig. 11
%                                                of the thesis
%
%   parms  42 slots (column). Slot 43 (f_l) is optional and absent: the core uses f_l = 0.01.
%   y0     31 state variables in the LABORATORY convention (compression positive).
%   TESTS  initial state of the four undrained cyclic tests Tx06-Tx09 (Saldana Sotelo 2023, MSc thesis,
%          Universidad de Concepcion, Table A.4.1; Saldana et al. 2024): e0 and Dr as in the laboratory
%          index used in the calibration (Table A.4.1 prints 0.701 and 0.671 for Tx06 and Tx07 and Dr
%          79.2 % for Tx07), p0 = first measured mean effective stress, CSR = q_amp/p0 as in the thesis
%          (twice the ratio q_amp/(2 p0) of Table A.4.1), N_test from Table A.4.1. TESTS.TMD holds the
%          initial state of the four drained tests of Fig. 11.
%
% The experimental data of Bio-Bio sand are NOT redistributed (see data/BBS/README.md). Everything
% else runs without them: the published tests start from the states tabulated here.
%
% parms(34) is ETA (rate of the cyclic memory) and parms(35) is KAPPA (rate of the monotonic memory);
% model/M_istr.m reads them as eta_rate_cyc and kappa_rate_mon. See docs/NOMENCLATURE.md.

if nargin < 1 || isempty(set), set = 'cyclic'; end
assert(any(strcmp(set, {'cyclic', 'monotonic'})), 'params_BBS: set must be ''cyclic'' or ''monotonic''.');

%% integration
max_ksub = 1000;   err_tol = 1.0e-4;   tol_NR = 1.0e-10;

%% hypoplasticity (von Wolffersdorff 1996) with the friction ramp
phi_c = 36.0;    h_s = 6.08e6;    n = 0.3020;
e_d0  = 0.596;   e_c0 = 0.958;    e_i0 = 1.150*0.958;      % = 1.1017
alpha = 0.48;    b_phi = 0.1;     r_KG = 1.13;   c = 0;

%% intergranular strain and cyclic memory
m_R = 1.62;   m_T = 1.10;   R = 1.0e-4;   beta_r = 0.15;   chi_0 = 5.5;   chi_max = 10;

%% fabric and semifluidized state
c_z = 300;   z_max = 12.0;
c_l = 40;    n_l = 8;    p_th = 10.0;    c_r = 1000;    lambda_1 = 0.07;    lambda_2 = 3.0;

%% memories
eta = 0.35;   kappa = 15;

%% switches and retired slots
Omega = 0;   a = 0;   p_inr = 0;   p_in = 0;   % a must stay 0 in this release (see params_KFS.m)
y_h = 1;   c_h = 0;   h_0 = 0.25;   g = 0;   k_chi = 0;   psi_ref = -0.218;   tent = 2;
q_br = 2.0;

z0 = zeros(6, 1);                              % the 11 m layer starts without fabric, z(0) = 0
if strcmp(set, 'monotonic')                    % drained tests of Cuitino et al. (2015) and Brown et al. (2024), thesis Fig. 11
    e_d0 = 0.610;   e_c0 = 0.990;   e_i0 = 1.1385;
    m_R = 2.2;   m_T = 1.1;   z_max = 8.0;   lambda_1 = 0.15;   eta = 20;
    z0 = 0.5*[-1/sqrt(6); -1/sqrt(6); 2/sqrt(6); 0; 0; 0];   % ||z(0)|| = 0.5 along triaxial compression
end
e0_doc = e_i0 - .71*(e_i0 - e_d0);             % not read by the core (documentary)

assert(y_h == 1 && g == 0 && k_chi == 0, ['params_BBS: y_h = 1, g = 0 and k_chi = 0 give the cyclic ' ...
       'memory of the thesis and the manuscript (detector y_delta, docs/NOMENCLATURE.md).']);
assert(a == 0, 'params_BBS: a must be 0 (p_in at reversals is not implemented in model/updateModel_stress.m).');
parms = [max_ksub; err_tol; tol_NR; phi_c; h_s; n; e_d0; e_c0; e_i0; alpha; b_phi; ...   %  1-11
         m_R; m_T; R; beta_r; chi_0; c; r_KG; q_br; ...                                    % 12-19
         c_z; z_max; c_l; n_l; p_th; c_r; lambda_1; lambda_2; e0_doc; chi_max; Omega; ...  % 20-30
         a; p_inr; p_in; eta; kappa; y_h; c_h; h_0; g; k_chi; psi_ref; tent];              % 31-42

%% initial state (laboratory convention); stress and void ratio are set by init_state.m
eps0   = zeros(6, 1);
sig0   = [100; 100; 100; 0; 0; 0];
e0     = e0_doc;
delta0 = -0.8*R/sqrt(3)*[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];                % ||delta(0)|| = 0.8 R
y0     = [eps0; sig0; e0; delta0; z0; zeros(6, 1)];

%% published undrained cyclic tests (Saldana Sotelo 2023, Table A.4.1)
%   N_test = cycles to r_u = 95 % reported by Saldana; ncyc = ceiling of cycles for the integration,
%   ceil(15*N_test) + 15 as in the thesis run (the integration stops one cycle after liquefaction)
%   CSR = q_amp/p0 (Table A.4.1 prints q_amp/(2 p0) = 0.197, 0.206, 0.325, 0.169)
%            id   e0     Dr     CSR    p0       ncyc  N_test
T = [         6   0.700  71.0   0.394  102.914  180   11
              7   0.670  79.0   0.412  104.064  150    9
              8   0.709  69.0   0.650  103.645   45    2
              9   0.685  75.5   0.338  103.701  315   20 ];
TESTS = struct('id', T(:,1), 'e0', T(:,2), 'Dr', T(:,3), 'CSR', T(:,4), 'p0', T(:,5), ...
               'q_amp', T(:,4).*T(:,5), 'ncyc', T(:,6), 'N_test', T(:,7));

%% drained monotonic tests of another Bio-Bio sample (Cuitino et al. 2015; Brown et al. 2024; thesis Fig. 11)
%   isotropic consolidation to p0 = 100 kPa, then drained triaxial compression to eps_a = 30 %
TESTS.TMD = struct('id', (1:4)', 'p0', 100*ones(4,1), 'e0', [0.890; 0.738; 0.688; 0.644]);
end
