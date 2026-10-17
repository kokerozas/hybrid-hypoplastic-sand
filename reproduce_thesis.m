% REPRODUCE_THESIS - reruns the element tests behind the main numbers of the thesis and compares them
% with the reference values in reference/reference_values.csv (thesis run of 2026-09-23, MATLAB R2020a
% 9.8.0.1323502, with this same constitutive core).
%
%   KFS TCU   6 undrained cyclic tests: cycles to initial liquefaction (half-wave counting)
%   BBS TCU   4 undrained cyclic tests: cycles to initial liquefaction (counting inside the cycle)
%   KFS TMD  25 drained monotonic tests: peak and final q, minimum volumetric strain
%   KFS TMU  14 undrained monotonic runs (12 distinct tests): minimum p', final q
%
% Not covered: the oedometer calibration of KFS, the drained tests of BBS (thesis Fig. 11) and the
% tests with constant strain amplitude.
%
% Each check is EXACT (relative difference <= 1e-12), MATCH (cycles to liquefaction within 0.01
% cycles, monotonic quantities within 1e-3 relative) or FAIL. On the platform of the reference run all
% checks are expected to be EXACT; on another release or processor a difference of one integration
% step is possible and is reported as MATCH. The running time is printed at the end.

clear variables; close all; clc
root = fileparts(mfilename('fullpath'));
assert(~isempty(root), 'reproduce_thesis: run the saved file (F5 or run(''reproduce_thesis''))');
addpath(root, fullfile(root, 'model'), fullfile(root, 'params'), fullfile(root, 'tools'));
for fn = {'input_data', 'init_state', 'hyp_stif', 'updateModel', 'updateModel_stress', ...
          'params_KFS', 'params_BBS', 'cyclic_test'}
    assert(strncmpi(which(fn{1}), root, numel(root)), ...
           'reproduce_thesis: %s.m is taken from %s (current folder or earlier on the path).', fn{1}, which(fn{1}));
end
global log_internal_vars log_index %#ok<GVMIS>
log_internal_vars = struct();   log_index = [];

RUN = {'KFS', 'TCU';  'BBS', 'TCU';  'KFS', 'TMD';  'KFS', 'TMU';  'KFS', 'TMU2'};   % delete rows to skip
TOL_EXACT = 1e-12;   TOL_N = 0.01;   TOL_REL = 1e-3;

fid = fopen(fullfile(root, 'reference', 'reference_values.csv'), 'r');
C = textscan(fid, '%s%s%f%s%f%s%q', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);
assert(numel(C{5}) == 113 && all(isfinite(C{5})) && all(isfinite(C{3})), ...
       'reproduce_thesis: reference/reference_values.csv was not read correctly.');
REF = struct('sand', {C{1}}, 'test', {C{2}}, 'id', C{3}, 'quantity', {C{4}}, 'value', C{5});

REPRODUCING = true;   % input_data.m keeps SAND, TEST and ID set here
BBS_SET = 'cyclic';   SIG1 = 15;   DX_C = 0.04;   DX_U = 0.01;
OUT = cell(0, 6);
t00 = tic;
try
    for r = 1:size(RUN, 1)
        SAND = RUN{r,1};   TEST = RUN{r,2};
        sel = strcmp(REF.sand, SAND) & strcmp(REF.test, TEST);
        ids = unique(REF.id(sel), 'stable')';
        fprintf('\n%s %s: %d tests\n', SAND, TEST, numel(ids));
        for ID = ids
            input_data;   init_state;
            t0 = tic;
            if strcmp(TEST, 'TCU')
                if strcmp(SAND, 'KFS'), NPOST = Inf; end       % the whole cycle budget, as in the thesis run
                [SS,EE,INV_S,INV_E,HARD,NLIQ] = cyclic_test(y0, parms, path_info, q_amp, NPOST);
                val = struct('N_liq_halfwave', NLIQ.halfwave, 'N_liq_cycle', NLIQ.cycle);
            else
                [SS,EE,INV_S,INV_E,HARD] = updateModel(y0, parms, nspb, path_info);
                q = SS(:,3) - SS(:,1);   p = INV_S(:,1);
                val = struct('q_peak', max(q), 'eps_v_min', min(INV_E(:,1)), 'q_end', q(end), 'p_min', min(p));
            end
            rows = find(sel & REF.id == ID)';
            msg = sprintf('  %s %-4s %3d  %5.1f s', SAND, TEST, ID, toc(t0));
            for k = rows
                qty = REF.quantity{k};   got = val.(qty);   ref = REF.value(k);
                if ref == 0, d = abs(got); else, d = abs(got - ref)/abs(ref); end   % absolute when ref = 0
                if strncmp(qty, 'N_liq', 5), match = abs(got - ref) <= TOL_N; else, match = d <= TOL_REL; end
                if d <= TOL_EXACT, st = 'EXACT'; elseif match, st = 'MATCH'; else, st = 'FAIL'; end
                OUT(end+1, :) = {SAND, TEST, ID, qty, d, st}; %#ok<SAGROW>
                msg = [msg sprintf('  | %s %.6g (ref %.6g) %s', qty, got, ref, st)]; %#ok<AGROW>
            end
            fprintf('%s\n', msg);
        end
    end
catch err
    clear REPRODUCING
    rethrow(err)
end
clear REPRODUCING

n = size(OUT, 1);   ne = nnz(strcmp(OUT(:,6), 'EXACT'));   nm = nnz(strcmp(OUT(:,6), 'MATCH'));
nf = n - ne - nm;
dall = [OUT{:,5}];   dmax = max(dall);   if any(isnan(dall)), dmax = NaN; end
fprintf('\n%d checks: %d EXACT, %d MATCH, %d FAIL; largest difference %.1e  (%.1f min)\n', ...
        n, ne, nm, nf, dmax, toc(t00)/60);
if nf == 0, fprintf('RESULT: PASS\n'), else, fprintf(2, 'RESULT: FAIL\n'); disp(OUT(strcmp(OUT(:,6), 'FAIL'), 1:5)), end
