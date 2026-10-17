function [SS,EE,INV_S,INV_E,HARD] = updateModel_stress(y0,parms,nspb,path_info,sig_input,sig_input2)

% UPDATE: compute the evolution of strain, stress and internal variables
%         over the prescribed loading path
%
% version 1.0 for rate-independent elastoplastic models
%
% integration methods implemented:
% explicit adaptive Runge-Kutta (3rd order) with error control 
%
% written by C. Tamagnini - jan. 2012

%% initialize load rate vector

V=zeros(6,1);

%% recover input data
code = [];
tol_f=parms(1);

ny = max(size(y0));

eps  = y0(1:6,1);
sig  = y0(7:12,1);
qint = y0(13:ny,1);

icode = path_info(1,:);
nstep = path_info(2,:);
DX    = path_info(3,:);


dx = DX./nstep;  % increment size for each sp branch
% dx = DX;  % increment size for each sp branch

%% copy initial state into global output arrays

kstep = 1;

[p,q,z] = inv_s(0 - sig);        % 2026-09-19 (rep. 152): inv_s/inv_e reciben tension y deformacion del MODELO (compresion negativa); p, q, eps_v, eps_s salen ingenieriles
[epsv,epss] = inv_e(0 - eps);

SS(kstep,:)=sig';
EE(kstep,:)=eps';
INV_S(kstep,:)=[p,q,z];
INV_E(kstep,:)=[epsv,epss];
HARD(kstep,:)=[qint'];

%% loop over the SPB

for i=1:nspb
    if -(-1)^i>0
    sig_obj = -sig_input*(-1)^i;
    else
    sig_obj = -sig_input2*(-1)^i;   
    end
    %disp(' ')
	%disp(['STRESS PATH BRANCH # ', int2str(i)])
	%disp(' ')

	n=100000; kstep0 = kstep;   % n: tope pasos/rama (SMOKE=6000; VEREDICTO=100000). kstep0: paso inicial de la rama
	k=icode(i);
    dload = dx(i);

% loop over load steps

% 	while abs(sig(6)-sig_obj) > 1e-3
	while abs(sig(3) - sig(1)-sig_obj) > 1e-3 && (kstep - kstep0) < n

		kstep=kstep+1;
		%disp(['STEP # ',int2str(kstep-1)])
 
% load vector V
        txx = min(max(abs(sig(3) - sig(1)-sig_obj) , 1),5);
%         txx = min(max(abs(sig(6) -sig_obj) , 1),5);
%         V(6)=dload*txx*((0.001-1)/100*p + 1);
        V(6)=txx*dload;
%         V(6)=txx*dload*exp(qint(end-4));
%         V(6)= dload*exp(qint(end-1));

% res   = sig_obj - (sig(3)-sig(1));
% scale = min(max(abs(res)/q_ref, vmin), vmax);  % e.g., q_ref=5 kPa, vmin=0.5, vmax=5
% V(6)  = sign(res) * dload * scale;


% initial state

    y_k=[eps',sig',qint']';
%  hypoplastic evolution equation

% ----- frontera laboratorio/modelo (2026-09-19, rep. 152): el nucleo integra en compresion
%       NEGATIVA (von Wolffersdorff 1996); y_k, V, SS, EE, HARD siguen en compresion POSITIVA.
%       espejo_estado es una involucion exacta (0 - x en las ranuras impares); 0 - V idem.
        [y_k , code] = hypo_update(espejo_estado(y_k),0 - V,parms,k);
        y_k = espejo_estado(y_k);
    
        if code == 1
            return
        end
% recover state variables and invariants at the end of substepping

		eps  = y_k(1:6,1);
		sig  = y_k(7:12,1);
		qint = y_k(13:ny,1);

        [p,q,z] = inv_s(0 - sig);        % 2026-09-19 (rep. 152): inv_s/inv_e reciben tension y deformacion del MODELO (compresion negativa); p, q, eps_v, eps_s salen ingenieriles
        [epsv,epss] = inv_e(0 - eps);
        
%         figure(3)
%         plot(p,sig(6),'ro') ;hold on; pause(0.0001)
% move final states into stress and strain path vectors

        SS(kstep,:)=sig';
        EE(kstep,:)=eps';
        INV_S(kstep,:)=[p,q,z];
        INV_E(kstep,:)=[epsv,epss];
        HARD(kstep,:)=[qint'];

        
%    simpleshear     
%         if -(-1)^i > 0
%           if sig_obj < sig(6) 
%               sig_obj - sig(6);
%               break
%           end
%         else
%             
%             if sig_obj > sig(6)
%                 sig_obj - sig(6);
%                 break
%             end
%         end
        
% % triaxial

        if -(-1)^i > 0
          if sig_obj < sig(3) - sig(1)
              sig(3) - sig(1) - sig_obj;
              break
          end
        else
            
            if sig_obj > sig(3) - sig(1)
                sig_obj - (sig(3) - sig(1));
                break
            end
        end
                
%         if abs(sig(3)-sig(1) - sig_obj) < 1
%             return
%         end
% bottom of the loop over load steps

	end

% bottom of the loop over stress path branches
 
end
