function [E,S] = constraints(k,y)

% CONSTRAINTS: evaluate constraint matrices E and S from
%              code k and current state (stored in y)
%              note: new cases can be added at the end of
%              the block-if

%% recover material state (in case E and S depend on sig and q)

ny = max(size(y));

sig = y(7:12,1);
q   = y(13:ny,1);

%% initialize E and S
 
E = zeros(6);
S = zeros(6);

%% build matrices E and S

if k==1

% strain controlled undrained TX compression (axis direction: x_3)

	E(1,1)=1;
	E(1,3)=0.5;
	E(2,2)=1;
	E(2,3)=0.5;
	E(3,4)=1;
	E(4,5)=1;
	E(5,6)=1;
    E(6,3)=1;

elseif k==2

% strain controlled ED compression (axis direction: x_3)

	E(1,1)=1;
	E(2,2)=1;
	E(3,4)=1;
	E(4,5)=1;
	E(5,6)=1;
        E(6,3)=1;
    
elseif k==3

% stress controlled drained TX compression (axis direction: x_3)

	S(1,1)=1;
	S(1,3)=0;
	S(2,2)=1;
	S(2,3)=0;
	S(3,4)=1;
	S(4,5)=1;
	S(5,6)=1;
	S(6,3)=1;

elseif k==4

% mixed control drained TX compression

	S(1,1)=1;
	S(2,2)=1;
	S(3,4)=1;
	S(4,5)=1;
	S(5,6)=1;
	E(6,3)=1;

elseif k==5

% mixed control ED compression

	E(1,1)=1;
	E(2,2)=1;
	E(3,4)=1;
	E(4,5)=1;
	E(5,6)=1;
	S(6,3)=1;

elseif k==6

% mixed control plane strain compression

	S(1,1)=1;
	S(3,4)=1;
	S(4,5)=1;
	S(5,6)=1;
	E(2,2)=1;
	E(6,3)=1;

elseif k==7

% strain controlled undrained simple shear

	E=eye(6);

elseif k == 8

% stress (mixed) controlled simple shear

   	S(3,3)=1;
	S(4,4)=1;
	S(5,5)=1;
	S(6,6)=1;
	E(1,1)=1;
	E(2,2)=1;
else
	E=eye(6);
    E(6,6) = 0;
    S(6,6) = 1;
end