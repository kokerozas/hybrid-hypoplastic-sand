% INIT_STATE: initialize the state variables
%
% Same role as init_state.m in the original driver: it builds the state vector y0 (31 x 1) in the
% LABORATORY convention (compression positive). The drivers in model/ map it to the convention of
% the constitutive equations (compression negative, von Wolffersdorff 1996) with model/espejo_estado.m.
%
%   y0(1:6)    strain eps            (Voigt order 11 22 33 12 23 13; shear strains engineering)
%   y0(7:12)   effective stress sig  [kPa]
%   y0(13)     void ratio e
%   y0(14:19)  intergranular strain delta (Niemunis & Herle 1997)
%   y0(20:25)  fabric tensor z (Liao et al. 2022)
%   y0(26)     state variable l of the semifluidized state
%   y0(28)     cyclic strain memory eps_cyc;   y0(30) monotonic strain memory eps_mon
%   y0(27), y0(29), y0(31) are not used by the equations
%
% delta(0) and z(0) are properties of each sand and come from params_KFS.m / params_BBS.m.

y0 = y0_sand;
y0(7:9) = p0;      % isotropic initial effective stress
y0(13)  = e0;      % initial void ratio
