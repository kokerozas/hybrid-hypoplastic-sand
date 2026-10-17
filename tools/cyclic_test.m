function [SS, EE, INV_S, INV_E, HARD, NLIQ] = cyclic_test(y0, parms, path_info, q_amp, NPOST)
%CYCLIC_TEST  undrained cyclic triaxial test with constant stress amplitude, cycle by cycle.
%
%   [SS, EE, INV_S, INV_E, HARD, NLIQ] = cyclic_test(y0, parms, path_info, q_amp, NPOST)
%
% Integrates the loading path of input_data.m two branches (one cycle) at a time with
% model/updateModel_stress.m and stops NPOST cycles after initial liquefaction (NPOST = Inf runs all
% the branches of path_info; NPOST >= 1 so that the half-wave count closes its half cycle). Cycle by
% cycle and in a single call the trajectory is the same by construction, because the complete state
% is passed between cycles (measured to the last bit for KFS TCUI10 and BBS Tx07); the loop is only
% needed to know where each cycle starts.
%
% Initial liquefaction is the first state with p'/p'_0 <= 0.05 (r_u >= 0.95); the half-wave count uses
% the strict form r_u > 0.95 of the thesis runs, which differs only at exact equality. The thesis counts
% the cycles in two ways, which differ only in the definition of the cycle axis N:
%   NLIQ.halfwave  N advances 0.5 at every reversal of q (tools/cycles_halfwave.m): KFS values
%   NLIQ.cycle     N = position of the state inside its cycle (the counter of the thesis for BBS)
% Both return the cycle of the same state; they differ by a fraction of a cycle.
%
% Outputs SS, EE, INV_S, INV_E, HARD as in updateModel.m (laboratory convention), one row per state,
% without the repeated first row of each cycle. NLIQ also holds both cycle axes (N_halfwave, N_cycle).

if nargin < 5 || isempty(NPOST), NPOST = 2; end
NPOST = max(NPOST, 1);
p0 = y0(7);
ncyc = floor(size(path_info, 2)/2);
y = y0;
C = cell(ncyc, 6);
Ncyc = NaN;
nc = 0;
for c = 1:ncyc
    nc = c;
    [S, E, I, IE, H] = updateModel_stress(y, parms, 2, path_info(:, 2*c-1:2*c), q_amp, q_amp);
    y = [E(end,:), S(end,:), H(end,:)]';
    n = size(S, 1);
    if isnan(Ncyc)
        ru = 1 - I(:,1)/p0;
        j = find(ru >= 0.95, 1, 'first');
        if ~isempty(j), Ncyc = c - 1 + (j - 1)/(n - 1); end
    end
    if c == 1, s = 1:n; else, s = 2:n; end           % row 1 repeats the end of the previous cycle
    C(c,:) = {S(s,:), E(s,:), I(s,:), IE(s,:), H(s,:), c - 1 + (s(:) - 1)/(n - 1)};
    qc = S(:,3) - S(:,1);                             % the cycle must reach +q_amp and then -q_amp; if not,
    if max(qc) < q_amp - 1e-3 || qc(end) > -q_amp + 1e-3   % hypo_update gave up or a branch hit the cap of
        warning('cyclic_test: cycle %d did not reach +-q_amp; integration ends.', c);   % updateModel_stress
        break
    end
    if ~isnan(Ncyc) && c - Ncyc >= NPOST, break, end
end
C = C(1:nc, :);
SS = vertcat(C{:,1});   EE = vertcat(C{:,2});   INV_S = vertcat(C{:,3});
INV_E = vertcat(C{:,4});   HARD = vertcat(C{:,5});

p = INV_S(:,1);   q = SS(:,3) - SS(:,1);
Nh = cycles_halfwave(q);
jh = find(1 - p/p(1) > 0.95, 1);                  % strict '>' as in the KFS counter of the thesis
NLIQ = struct('halfwave', NaN, 'cycle', Ncyc, 'N_halfwave', Nh, 'N_cycle', vertcat(C{:,6}), 'cycles_run', nc);
if ~isempty(jh), NLIQ.halfwave = Nh(jh); end
end
