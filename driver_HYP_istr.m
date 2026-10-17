% DRIVER_HYP_ISTR - element-test driver of the hybrid hypoplastic model for sands
%
% Constitutive driver for the numerical integration of hypoplastic constitutive equations with an
% explicit adaptive Runge-Kutta scheme and the mixed-control algorithm of Bardet & Choucair (1991),
% following the driver written by C. Tamagnini (2012) with the modifications for hypoplasticity by
% D. Masin (course material of the ALERT Olek Zienkiewicz School 2012, Dresden). The workflow is the
% original one:
%
%   original header: written by C. Tamagnini - Jan. 2012; modifications for hypoplasticity by D. Masin
%
%   input_data   material constants (parms) and loading path (nspb, path_info)
%   init_state   initial state vector (y0)
%   updateModel  integration over the loading path (updateModel_stress for stress amplitudes)
%   xyplot_vmh   plots
%
% The constitutive equations in model/ are those of the thesis: hypoplasticity of von Wolffersdorff
% (1996), intergranular strain of Niemunis & Herle (1997), fabric of Liao et al. (2022),
% semifluidized state of Barrero et al. (2020) in the form of Liao et al. (2022), and the cyclic and
% monotonic strain memories with the friction ramp.
%
% Usage: edit block 1 of input_data.m (sand, test, published test or your own state) and run this
% script. To rerun the main numbers of the thesis (113 values, see README.md), run reproduce_thesis.m.

clear variables; close all; clc
root = fileparts(mfilename('fullpath'));
assert(~isempty(root), 'driver_HYP_istr: run the saved file (F5 or run(''driver_HYP_istr''))');
addpath(root, fullfile(root, 'model'), fullfile(root, 'params'), fullfile(root, 'tools'));
for fn = {'input_data', 'init_state', 'xyplot_vmh', 'hyp_stif', 'updateModel', 'updateModel_stress', ...
          'params_KFS', 'params_BBS', 'cyclic_test'}
    assert(strncmpi(which(fn{1}), root, numel(root)), ...
           'driver_HYP_istr: %s.m is taken from %s (current folder or earlier on the path).', fn{1}, which(fn{1}));
end
global log_internal_vars log_index %#ok<GVMIS>
log_internal_vars = struct();   log_index = [];   % internal logging of the core: off

disp('===========================================')
disp('                DRIVER HYP                 ')
disp('   hybrid hypoplastic model for sands      ')
disp('===========================================')
disp('  driver: C. Tamagnini (2012), with the    ')
disp('  modifications for hypoplasticity by      ')
disp('  D. Masin (ALERT School 2012, Dresden)    ')
disp('===========================================')

%% input data and initial material state
input_data;
init_state;

%% call integration procedure
fprintf('\n%s %s %s   p0 %.3f kPa   e0 %.4f\n', SAND, TEST, num2str(ID), p0, e0);
t0 = tic;
if strcmp(TEST, 'TCU')
    fprintf('q_amp %.4f kPa, up to %d cycles ...\n', q_amp, nspb/2);
    [SS,EE,INV_S,INV_E,HARD,NLIQ] = cyclic_test(y0, parms, path_info, q_amp, NPOST);
    Nliq = NLIQ.(NLIQ_RULE);
    if isnan(Nliq)
        fprintf('no initial liquefaction (p''/p''0 <= 0.05) in %d cycles', NLIQ.cycles_run);
    else
        fprintf('N_liq = %.3f cycles (first p''/p''0 <= 0.05; counting: %s)', Nliq, NLIQ_RULE);
    end
    if ~isempty(ID), fprintf('   test: %g', TESTS.N_test(TESTS.id == ID)); end
    fprintf('\n');
else
    [SS,EE,INV_S,INV_E,HARD] = updateModel(y0, parms, nspb, path_info);
    if size(SS, 1) < 1 + sum(path_info(2,:))
        warning('driver_HYP_istr: updateModel stopped at step %d of %d (substepping failure).', ...
                size(SS, 1) - 1, sum(path_info(2,:)));
    end
    q = SS(:,3) - SS(:,1);
    fprintf('final state: eps_a %.2f %%   q %.2f kPa   p %.2f kPa   e %.4f\n', ...
            100*EE(end,3), q(end), INV_S(end,1), HARD(end,1));
end
fprintf('%.1f s\n', toc(t0));

%% output section: plots (and, if wanted, the arrays as text files, as in the original driver)
xyplot_vmh;
% save inv_s.txt INV_S -ascii;   save inv_e.txt INV_E -ascii;   save stresses.txt SS -ascii;
% save strains.txt EE -ascii;    save hard.txt HARD -ascii;
