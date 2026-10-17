% INPUT_DATA: choose the test, read the material constants and set the loading path
%
% Same role as input_data.m in the original constitutive driver (C. Tamagnini 2012, with the
% modifications for hypoplasticity of D. Masin): it leaves in the workspace the parameter vector
% parms and the loading path
%
%   nspb          = no. of stress path branches
%   path_info     = [icode; nstep; DX], one column per branch
%   icode(1:nspb) = path type code of each branch (model/constraints.m)
%   nstep(1:nspb) = number of steps of each branch
%   DX(1:nspb)    = total "load" increment of each branch, the quantity in row 6 of the constraints:
%                   axial strain eps_33 (absolute, not %) for codes 1, 2, 4 and 6, axial stress
%                   sigma_33 [kPa] for codes 3 and 5, shear strain gamma_13 for code 7 and shear
%                   stress sigma_13 [kPa] for code 8
%
% loading path codes (axis direction x_3)
%   1. strain controlled undrained TX       5. mixed control ED compression
%   2. strain controlled ED compression     6. mixed control plane strain compression
%   3. stress controlled drained TX         7. strain controlled undrained simple shear
%   4. mixed control drained TX             8. stress (mixed) controlled simple shear
%
% TCU (undrained cyclic): every branch is code 1 and ends when the deviator stress q reaches
% +q_amp or -q_amp (model/updateModel_stress.m); DX only sets the size of the strain steps.

%% 1. what to run  ------------------------------------------------------ edit this block
if ~exist('REPRODUCING', 'var')
    SAND = 'KFS';   % 'KFS' Karlsruhe fine sand | 'BBS' Bio-Bio sand
    TEST = 'TCU';   % 'TCU' undrained cyclic, constant stress amplitude
                    % 'TMU' undrained monotonic (KFS also 'TMU2', pressure series)
                    % 'TMD' drained monotonic | 'OED' oedometer, loading and unloading
    ID   = 10;      % a published test (list below) or [] for your own state and loading
    % your own test (ID = []):
    p0 = 100;   e0 = 0.80;              % initial isotropic effective stress [kPa] and void ratio
    q_amp = 30;   NCYC = 20;            % TCU: deviator stress amplitude [kPa] and number of cycles
    EPS_A = 0.30;                       % TMU, TMD: final axial strain (0.30 = 30 %; TMU < 0: extension)
    SIG1 = 15;   DX_C = 0.04;   DX_U = 0.01;   % OED: initial stress [kPa], axial strain loading and unloading
    BBS_SET = 'cyclic';                 % BBS: 'cyclic' (11 m layer) | 'monotonic' (drained tests, thesis Fig. 11)
end
% Published tests (ID):
%   KFS  TCU 8 9 10 18 19 20 | TMD 1-25 | TMU 1-6 | TMU2 1-8 | OED 1-12
%        data of Wichtmann & Triantafyllidis (2016), included in data/KFS (OED: loaded up to the
%        maximum measured stress; the thesis calibrated C_c and C_s with a per-specimen recipe)
%   BBS  TCU 6 7 8 9 (Tx06-Tx09, Saldana Sotelo 2023) | TMD 1-4 (Cuitino et al. 2015; Brown et al. 2024; thesis Fig. 11):
%        initial states in params/params_BBS.m; the measured curves are not redistributed (data/BBS)

%% 2. material constants  ------------------------------------------------------------------
root = fileparts(mfilename('fullpath'));
DATA = [];                       % measured curves of the published test, when available (not read by reproduce_thesis)
switch SAND
    case 'KFS'
        if strcmp(TEST, 'TCU'), SET = 'cyclic'; else, SET = 'monotonic'; end
        [parms, y0_sand, TESTS] = params_KFS(SET);
    case 'BBS'
        if strcmp(TEST, 'TCU'), SET = 'cyclic';
        elseif strcmp(TEST, 'TMD') && ~isempty(ID), SET = 'monotonic';   % the published drained tests
        else, SET = BBS_SET;
        end
        [parms, y0_sand, TESTS] = params_BBS(SET);
    otherwise
        error('input_data: SAND must be ''KFS'' or ''BBS''.')
end
dat = fullfile(root, 'data', SAND);

%% 3. initial state of the published test and loading path  ---------------------------------
NSTEP = 250;   DX1 = 0.0005;   DX = 0.001;   % TCU: strain step DX/NSTEP (first branch DX1/NSTEP); each branch runs to q = +-q_amp
NPOST = 2;                                   % TCU: cycles integrated after liquefaction (Inf = all NCYC); BBS: 1
NLIQ_RULE = 'halfwave';                      % TCU: how cycles are counted (tools/cyclic_test.m); BBS: 'cycle'
switch TEST
    case 'TCU'
        if ~isempty(ID)
            k = find(TESTS.id == ID, 1);
            assert(~isempty(k), 'input_data: %s TCU %d is not a published test.', SAND, ID);
            p0 = TESTS.p0(k);   e0 = TESTS.e0(k);   q_amp = TESTS.q_amp(k);   NCYC = TESTS.ncyc(k);
            if strcmp(SAND, 'KFS')
                f = fullfile(dat, 'TCU', sprintf('TCUI%d.dat', ID));
                if isfile(f) && ~exist('REPRODUCING', 'var'), A = importdata(f);   DATA = struct('p', A.data(:,7), 'q', A.data(:,8), 'e1', A.data(:,9)); end
            else                                       % optional, not distributed: data/BBS/README.md
                f = fullfile(dat, 'TCU - Hector Saldana', sprintf('CSR%d.dat', ID - 5));
                if isfile(f) && ~exist('REPRODUCING', 'var')
                    A = readmatrix(f, 'FileType', 'text');
                    if size(A, 2) < 3, A = readmatrix(f, 'FileType', 'text', 'Delimiter', ','); end
                    A = A(all(isfinite(A(:,1:3)), 2), :);
                    A = A(find(A(:,2) > 0, 1):end, :);   DATA = struct('e1', A(:,1), 'p', A(:,2), 'q', A(:,3));
                end
            end
        end
        if strcmp(SAND, 'BBS'), NLIQ_RULE = 'cycle'; NPOST = 1; end   % the counters of the thesis
        assert(NCYC >= 1 && NCYC == round(NCYC), 'input_data: NCYC must be a positive integer.');
        nspb = 2*NCYC;
        path_info = [ones(1, nspb); NSTEP*ones(1, nspb); [DX1, DX*(-1).^(1:nspb-1)]];

    case {'TMD', 'TMU', 'TMU2'}
        NSTEP_M = 500;                                % steps of the monotonic branch (thesis: KFS 500, BBS 600)
        if strcmp(SAND, 'BBS'), NSTEP_M = 600; end
        if ~isempty(ID) && strcmp(SAND, 'BBS')
            assert(strcmp(TEST, 'TMD') && ismember(ID, TESTS.TMD.id), 'input_data: the published monotonic tests of BBS are TMD 1-4.');
            p0 = TESTS.TMD.p0(ID);   e0 = TESTS.TMD.e0(ID);   EPS_A = 0.30;
            f = fullfile(dat, 'TMD_biobio_cuitino', sprintf('DR%d.dat', ID));   % optional, not distributed
            if isfile(f) && ~exist('REPRODUCING', 'var')
                A = readmatrix(f, 'FileType', 'text');   A = A(all(isfinite(A(:,1:4)), 2), :);
                DATA = struct('e1', A(:,1), 'q', A(:,2), 'ev', A(:,4), 'p', p0 + A(:,2)/3);
            end
        elseif ~isempty(ID)
            switch TEST
                case 'TMD',  d = 'TMD';    ix = 'TMD_index.txt';    tp = 1;
                case 'TMU',  d = 'TMU';    ix = 'TMU_index.txt';    tp = [1 1 1 -1 -1 -1];
                case 'TMU2', d = 'TMU2';   ix = 'TMU2_index.txt';   tp = [1 1 1 1 -1 -1 -1 -1];
            end
            IX = read_index(fullfile(dat, d, ix));          % [-, p0, e0], one row per test
            assert(ID >= 1 && ID <= size(IX, 1) && ID == round(ID), 'input_data: KFS %s %d is not a published test.', TEST, ID);
            p0 = IX(ID,2);   e0 = IX(ID,3);
            if strcmp(TEST, 'TMD'), EPS_A = 0.30; else, EPS_A = tp(ID)*0.12; end   % as in the thesis
            f = fullfile(dat, d, sprintf('%s%d.dat', TEST(1:3), ID));
            if isfile(f) && ~exist('REPRODUCING', 'var')
                A = importdata(f);
                if strcmp(TEST, 'TMD'), DATA = struct('e1', A.data(:,1), 'ev', A.data(:,2), 'q', A.data(:,6), 'p', A.data(:,7));
                else,                   DATA = struct('e1', A.data(:,1), 'p', A.data(:,7), 'q', A.data(:,8));
                end
            end
        end
        if strcmp(TEST, 'TMD'), icode = 4; else, icode = 1; end
        nspb = 1;
        path_info = [icode; NSTEP_M; EPS_A];

    case 'OED'
        assert(SIG1 > parms(24), 'input_data: OED starts above p_th = %g kPa (semifluidized state).', parms(24));
        if ~isempty(ID)
            assert(strcmp(SAND, 'KFS'), 'input_data: published oedometer tests are those of KFS.');
            f = fullfile(dat, 'OED', sprintf('OE%d.dat', ID));
            assert(isfile(f), 'input_data: KFS OED %d is not a published test.', ID);
            fid = fopen(f, 'r');   Cf = textscan(fid, '%f%f%f', 'HeaderLines', 2);   fclose(fid);
            DATA = struct('s1', Cf{1}, 'e', Cf{3});         % columns: sigma_1 [kPa], eps_1 [%], e
            jc = find(DATA.s1 == max(DATA.s1), 1);
            e0 = interp1(DATA.s1(1:jc), DATA.e(1:jc), SIG1);  % void ratio of the test at the initial stress
            assert(isfinite(e0), 'input_data: SIG1 is outside the loading branch of OE%d.', ID);
            % loading up to the maximum measured stress: axial strain where the model reaches it
            y_try = y0_sand;   y_try(7:9) = SIG1;   y_try(13) = e0;
            [S_try, E_try] = updateModel(y_try, parms, 1, [2; 500; DX_C]);
            jt = find(S_try(:,3) >= max(DATA.s1), 1);
            if ~isempty(jt), DX_C = E_try(jt,3); end
        end
        p0 = SIG1;
        nspb = 2;
        path_info = [2 2; 500 2000; DX_C -DX_U];

    otherwise
        error('input_data: TEST must be ''TCU'', ''TMU'', ''TMU2'', ''TMD'' or ''OED''.')
end
