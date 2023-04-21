function [M,Ovals] = TDVP_1site (M,Hs,O,Nkeep,dt)
% < Description >
%
% [M,Ovals] = TDVP_1site (M,Hs,O,Nkeep,dt)
%
% Time-dependent variational principle (TDVP) method for simulating
% real-time evolution of matrix product state (MPS). The expectation
% values of local operator O for individual sites are evaluated for
% discrete time instances. The 1D chain system is described by the matrix
% product operator (MPO) Hamiltonian Hs.
%
% < Input >
% M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
%       defines the chain length. The leg convention of M{n} is as follows:
%
%    1      3   1      3         1        3
%   ---M{1}---*---M{2}---* ... *---M{end}---
%       |          |                 |
%       ^2         ^2                ^2
%
% Hs : [cell] MPO description of the Hamiltonian. Each Hs{n} acts on site n
%       , and is rank-4 tensor. The order of legs of Hs{n} is left-bottom-
%       right-top, where bottom (top) leg is to be contracted with bra
%       (ket) tensor:
%
%       |4          |4
%    1  |   3    1  |   3
%   ---Hs{1}---*---Hs{2}---*--- ...
%       |           |
%       |2          |2
%
% O : [matrix] Rank-2 tensor as a local operator acting on a site. The
%       expectation value of this operator at each chain site is to be
%       computed; see the description of the output 'Ovals' for detail.
% Nkeep : [integer] Maximum bond dimension.
% dt : [numeric] Real time step size. Each real-time evolution by step dt
%       consists of one pair of sweeps (left-to-right and right-to-left).
%
% < Output >
% M : [cell] The final MPS after real-time evolution.
% Ovals : [matrix] Ovals(1,n) indicates the expectation value of local 
%       operator O (input) at the site n after the time step dt (input).
%
% Written by S.Lee (Jun.13,2019): Written for SoSe 2019.
% Updated by S.Lee (Jun.13,2019): Revised for SoSe 2020.
% Updated by J.Shim (Jun.25.2022): Revised for SoSe 2022.


% tobj = tic2;

N = numel(M);

% % sanity check for input
if N < 2
    error('ERR: chain is too short.');
elseif numel(M) ~= numel(Hs)
    error('ERR: M has different lengths from that of Hs.');
elseif ~ismatrix(O)
    error('ERR: local operator O should be rank 2.');
end

for itN = (1:N)
    if size(Hs{itN},2) ~= size(Hs{itN},4)
        error(['ERR: The second and fourth legs of Hs{', ...
            sprintf('%i',itN),'} have different dimensions.']);
    elseif size(Hs{itN},2) ~= size(M{itN},2)
        error(['ERR: The second legs of Hs{', ...
            sprintf('%i',itN),'} and M{',sprintf('%i',itN),'} have different dimensions.']);
    end
end
% % %

% results
Ovals = zeros(1,N);

% % Hamiltonian for the left/right parts of the chain
Hlr = cell(1,N+2);
% Hlr{1} and Hlr{end} are dummies; they will be kept empty. These dummies
% are introduced for convenience.

% Since M is in right-canonical form by now, Hlr{..} are the right parts of
% the Hamiltonian. That is, Hlr{n+1} is the right part of Hamiltonian which
% is obtained by contracting M(n:end) with Hs(n:end). (Note the index for
% Hlr is n+1, not n, since Hlr{1} is dummy.)
for itN = (N:-1:1)
    T = permute(M{itN},[3 2 1]); % permute left<->right, to make use of updateLeft
    if itN == N
        % "remove" the right leg (3rd leg) of Hs{itN}, which is dummy,
        % by permuting to the last; MATLAB automatically suppresses the
        % trailing singleton dimensions
        H2 = permute(Hs{itN},[2 1 4 3]); % bottom-left-top (-right)
        Hlr{itN+1} = updateLeft([],[],T,H2,3,T);
    elseif itN == 1
        Hlr{itN+1} = [];
    else
        % permute left<->right, to make use of updateLeft
        H2 = permute(Hs{itN},[3 2 1 4]); % right-bottom-left-top
        Hlr{itN+1} = updateLeft(Hlr{itN+2},3,T,H2,4,T);
    end
end

% left -> right
for itN = (1:(N-1))
    % time evolution of site-canonical tensor ("A tensor") M{itN}, via
    % TDVP_1site_expHA
    
    Anew = TDVP_1site_expHA (Hlr{itN},Hs{itN},Hlr{itN+2},M{itN},dt/2);
    
    % update M{itN} and generate Cold by using Anew, via SVD
    [M{itN},S2,V2] =  ...
        svdTr(Anew,3,[1 2],Nkeep,-1); % set Stol as -1, not to truncate even zero singular values
    Cold = contract(diag(S2),2,2,V2,2,1);
    
    % update Hlr{itN+1} in effective basis
    if itN == 1
        % "remove" the left leg (1st leg) of Hs{itN}, which is dummy,
        % by permuting to the last; MATLAB automatically suppresses the
        % trailing singleton dimensions
        H2 = permute(Hs{itN},[2 3 4 1]); % bottom-right-top (-left)
        Hlr{itN+1} = updateLeft([],[],M{itN},H2,3,M{itN});
    else
        Hlr{itN+1} = updateLeft(Hlr{itN},3,M{itN},Hs{itN},4,M{itN});
    end
    
    % inverse time evolution of C tensor (Cold -> Cnew), via TDVP_1site_expHC
    Cnew = TDVP_1site_expHC (Hlr{itN+1},Hlr{itN+2},Cold,dt/2);
    
    % absorb Cnew into M{itN+1}
    M{itN+1} = contract(Cnew,2,2,M{itN+1},3,1);
end

itN = N; % right end
M{itN} = TDVP_1site_expHA (Hlr{itN},Hs{itN},Hlr{itN+2},M{itN},dt);

% right -> left
for itN = ((N-1):-1:1)
    % update M{itN+1} and generate Cold via SVD
    [U2,S2,M{itN+1}] = ...
        svdTr(M{itN+1},3,1,Nkeep,-1); % set Stol as -1, not to truncate even zero singular values
    Cold = contract(U2,2,2,diag(S2),2,1);
    
    % update Hlr{itN+2} in effective basis
    T = permute(M{itN+1},[3 2 1]); % permute left<->right, to make use of updateLeft
    if (itN+1) == N
        % "remove" the right leg (3rd leg) of Hs{N}, which is dummy,
        % by permuting to the last; MATLAB automatically suppresses the
        % trailing singleton dimensions
        H2 = permute(Hs{itN+1},[2 1 4 3]); % bottom-left-top (-right)
        Hlr{itN+2} = updateLeft([],[],T,H2,3,T);
    else
        % permute left<->right, to make use of updateLeft
        H2 = permute(Hs{itN+1},[3 2 1 4]); % right-bottom-left-top
        Hlr{itN+2} = updateLeft(Hlr{itN+3},3,T,H2,4,T);
    end
    
    % inverse time evolution of C tensor (Cold -> Cnew), via TDVP_1site_expHC
    Cnew = TDVP_1site_expHC (Hlr{itN+1},Hlr{itN+2},Cold,dt/2);
    
    % absorb Cnew into M{itN}
    M{itN} = contract(M{itN},3,3,Cnew,2,1);
    
    % time evolution of site-canonical tensor ("A tensor") M{itN}, via TDVP_1site_expHA
    M{itN} = TDVP_1site_expHA (Hlr{itN},Hs{itN},Hlr{itN+2},M{itN},dt/2);
end

% Measurement of local operators O; currently M is in site-canonical
% with respect to site 1
MM = []; % contraction of bra/ket tensors from left
for itN = (1:N)
    Ovals(1,itN) = trace(updateLeft(MM,2,M{itN},O,2,M{itN}));
    MM = updateLeft(MM,2,M{itN},[],[],M{itN});
end

% toc2(tobj,'-v');
% chkmem;


end


function Anew = TDVP_1site_expHA (Hleft,Hloc,Hright,Aold,dt)
% Time evolution for "A tensor" Aold by time step dt, by using Hlr tensors
% that are the Hamiltonian in effective basis. Anew is the result after
% time evolution.
% This subfunction is adapted from the Lanczos routine 'DMRG/eigs_1site.m'.
% The difference from 'eigs_1site' is that it considers the time evolution,
% not the ground state, and does not consider the orthonormal state.

% default parameters
N = 5;
minH = 1e-10;

Asz = [size(Aold,1),size(Aold,2),size(Aold,3)]; % size of ket tensor
Akr = zeros(numel(Aold),N+1); % Krylov vectors (vectorized tensors)
Akr(:,1) = Aold(:)/norm(Aold(:)); % normalize Aold

% In the Krylov basis, the Hamiltonian becomes tridiagonal
ff = zeros(N,1); % 1st diagonal
gg = zeros(N+1,1); % main diagonal

for itN = (1:(N+1))
    % contract Hamiltonian with ket tensor
    Atmp = TDVP_1site_HA(Hleft,Hloc,Hright,reshape(Akr(:,itN),Asz));
    Atmp = Atmp(:); % vectorize
    
    gg(itN) = Akr(:,itN)'*Atmp; % diagonal element; "on-site energy"
    
    if itN < (N+1)
        % orthogonalize Atmp w.r.t. the previous ket tensors
        Atmp = Atmp - Akr(:,(1:itN))*(Akr(:,(1:itN))'*Atmp);
        Atmp = Atmp - Akr(:,(1:itN))*(Akr(:,(1:itN))'*Atmp); % twice, to reduce numerical noise
        
        % norm
        ff(itN) = norm(Atmp);
        
        if ff(itN) > minH
            Akr(:,itN+1) = Atmp/ff(itN);
        else
            % stop iteration; truncate ff, gg
            ff(itN:end) = [];
            gg(itN+1:end) = [];
            Akr(:,(itN+1):end) = [];
            break;
        end
    end
end

% Hamiltonian in the Krylov basis
Hkr = diag(ff,1);
Hkr = Hkr + Hkr' + diag(gg);
[Vkr,Ekr] = eig((Hkr+Hkr')/2);

Anew = Akr*(Vkr*(diag(exp((-1i*dt)*diag(Ekr)))*Vkr(1,:).'));
Anew = Anew/norm(Anew); % normalize
Anew = reshape(Anew,Asz); % reshape back to rank-3 tensor

end


function Aout = TDVP_1site_HA (Hleft,Hloc,Hright,Ain)
% Adapted from the subfunction 'eigs_1site_HA' in 'DMRG/eigs_1site.m'.

% set empty tensors as 1, for convenience
if isempty(Hleft)
    Hleft = 1;
end
if isempty(Hright)
    Hright = 1;
end

Aout = contract(Hleft,3,3,Ain,3,1);
Aout = contract(Aout,4,[2 3],Hloc,4,[1 4]);
Aout = contract(Aout,4,[2 4],Hright,3,[3 2]);
    
end



function Anew = TDVP_1site_expHC (Hleft,Hright,Aold,dt)
% Time evolution for "C tensor" Aold by time step -dt, by using Hlr tensors
% that are the Hamiltonian in effective basis. Cnew is the result after
% time evolution.
% This subfunction is adapted from the subfunction 'TDVP_1site_expHA'
% above. The differences here are that the ket tensor is rank-2; that there
% is no the local Hamiltonian 'Hloc'; and that the time evolution is in an
% inverse direction.

% default parameters
N = 5;
minH = 1e-10;

Asz = [size(Aold,1),size(Aold,2)]; % size of ket tensor
Akr = zeros(numel(Aold),N+1); % Krylov vectors (vectorized tensors)
Akr(:,1) = Aold(:)/norm(Aold(:)); % normalize Aold

% In the Krylov basis, the Hamiltonian becomes tridiagonal
ff = zeros(N,1); % 1st diagonal
gg = zeros(N+1,1); % main diagonal

for itN = (1:(N+1))
    
    % contract Hamiltonian with ket tensor
    Atmp = TDVP_1site_HC(Hleft,Hright,reshape(Akr(:,itN),Asz));
    Atmp = Atmp(:); % vectorize
    
    gg(itN) = Akr(:,itN)'*Atmp; % diagonal element; "on-site energy"
    
    if itN < (N+1)
        % orthogonalize Atmp w.r.t. the previous ket tensors
        Atmp = Atmp - Akr(:,(1:itN))*(Akr(:,(1:itN))'*Atmp);
        Atmp = Atmp - Akr(:,(1:itN))*(Akr(:,(1:itN))'*Atmp); % twice, to reduce numerical noise
        
        % norm
        ff(itN) = norm(Atmp);
        
        if ff(itN) > minH
            Akr(:,itN+1) = Atmp/ff(itN);
        else
            % stop iteration; truncate ff, gg
            ff(itN:end) = [];
            gg(itN+1:end) = [];
            Akr(:,(itN+1):end) = [];
            break;
        end
    end
end

% Hamiltonian in the Krylov basis
Hkr = diag(ff,1);
Hkr = Hkr + Hkr' + diag(gg);
[Vkr,Ekr] = eig((Hkr+Hkr')/2);

Anew = Akr*(Vkr*(diag(exp((+1i*dt)*diag(Ekr)))*Vkr(1,:).')); % inverse
Anew = Anew/norm(Anew); % normalize
Anew = reshape(Anew,Asz); % reshape back to rank-2 tensor

end


function Aout = TDVP_1site_HC (Hleft,Hright,Ain)
% Adapted from the subfunction 'TDVP_1site_HA' above.

% set empty tensors as 1, for convenience
if isempty(Hleft)
    Hleft = 1;
end
if isempty(Hright)
    Hright = 1;
end

Aout = contract(Hleft,3,3,Ain,2,1);
Aout = contract(Aout,3,[2 3],Hright,3,[2 3]);
    
end
