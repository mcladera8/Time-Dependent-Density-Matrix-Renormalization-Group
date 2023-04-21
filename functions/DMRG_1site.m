function [M,E0,Eiter] = DMRG_1site(Hs,Minit,Nkeep,Nsweep,varargin)
% < Description >
%
% [M,E0,Eiter] = DMRG_1site(Hs,Minit,Nkeep,Nsweep [, option])
%
% Single-site density-matrix renormalization group (DMRG) 
% calculation to search for the ground state and its energy
% of a one-dimensional system, whose Hamiltonian is given 
% by the matrix product operator Hs.
%
% < Input >
% Hs : [1 x N cell array] Matrix product operator (MPO) of the
%	Hamiltonian. Each Hs{n} acts on site n, and is a rank-4
%	tensor. The order of legs of Hs{n} is 
%	left-bottom-right-top, where bottom (top) leg is to be 
%	contracted with bra (ket) tensor. The length N of Hs, 
%	i.e., numel(Hs), determines the chain length.
% Minit : [1 x N cell array] Initial MPS from which to start the 
%	ground state search
% Nkeep : [numeric] Maximum bond dimension of the matrix product 
%	state (MPS) to consider.
% Nsweep : [numeric] Number of sweeps will be 2*Nsweep, as there 
%	are Nsweep times of round trip 
%	(right -> left, left -> right).
%
% < Option >
% 'Econv',.. : [numeric] Convergence criterion for energy. If 
%	Einit - Efin < Econv, stop sweeping even if less than
%	Nsweep sweeps have been done so far. Here, Einit and
%	Efin are the energies before and after one 
%	(right -> left, left -> right) round trip, respectively.
%	(Default: -inf, i.e. no energy convergence criterion.)
%
% < Output >
% M : [1 x N cell array] The result MPS which is obtained 
%	variationally to have the minimum expectation value of the
%	Hamiltonian H. It is in *left-canonical* form, since the 
%	last sweep is from left to right.
% E0 : [numeric] The energy of M.
% Eiter : [N x (2*Nsweep) numeric array] Each element Eiter(m,n) 
%	means the variational energy in the m-th iteration within 
%	the n-th sweep. Odd n is for right-to-left sweep and even n 
%	for left-to-right sweep. Note that the iteration index m 
%	matches with the site index for the left-to-right sweep; 
%	the iteration m corresponds to the site (N+1-m) for 
%	the right-to-left sweep.
%
% Written by S.Lee (May 28,2019)
% Updated by S.Lee (May 23,2020): Revised for SoSe2020.
% Update by J.Shim (May 25.2022): Revised for SoSe2022.
% Update by J.Shim (Jun 07.2022): Introducing option 'Econv'.

tobj = tic2;

% % sanity check for input and option
N = numel(Hs);

if N < 2
    error('ERR: chain is too short.');
end

for itN = (1:N)
    if size(Hs{itN},2) ~= size(Hs{itN},4)
        error(['ERR: The second and fourth legs of Hs{', ...
            sprintf('%i',itN),'} have different dimensions.']);
    end
end
% % %

% default parameter
Econv = -inf;

% parsing option
while ~isempty(varargin)
    switch varargin{1}
        case 'Econv'
            Econv = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: check input!');
    end
end

% show message
fprintf('Single-site DMRG: search for the ground state\n');
fprintf(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', # of sweeps = ',sprintf('%i',Nsweep),' x 2\n']);

% initialize MPS by Minit
M = Minit;
M = canonForm(M,0); % bring into right-canonical form
M = canonForm(M,N); % bring into left-canonical form

% ground-state energy for each iteration
Eiter = nan(N,2*Nsweep);
% later, Eiter(end,end) will be taken as the final result E0

% % Hamiltonian for the left/right parts of the chain
Hlr = cell(1,N+2);


% Since M is in left-canonical form by now, Hlr{..} are the left parts of
% the Hamiltonian. That is, Hlr{n+1} is the left part of Hamiltonian which
% is obtained by contracting M(1:n) with Hs(1:n). (Note the index for Hlr
% is n+1, not n, since Hlr{1} is dummy.)
for itN = (1:N)
    if itN == 1
        % "remove" the left leg (1st leg) which is dummy by permuting to the last
        H2 = permute(Hs{itN},[2 3 4 1]);
        Hlr{itN+1} = updateLeft([],[],M{itN},H2,3,M{itN});
    else
        Hlr{itN+1} = updateLeft(Hlr{itN},3,M{itN},Hs{itN},4,M{itN});
    end
end

for itS = (1:Nsweep)
    % right -> left
    for itN = (N:-1:1)
        % Use eigs_1site to obtain the variationally chosen ket tensor
        % Ceff and energy expectation value Eeff
        [Ceff,Eeff] = eigs_1site(Hlr{itN},Hs{itN},Hlr{itN+2},M{itN});
        
        Eiter(N+1-itN,2*itS-1) = Eeff;
        
        % update M{itN} and M{itN-1} by using Ceff, via SVD
        % decompose Ceff
        [UT,ST,M{itN}] = svdTr(Ceff,3,1,Nkeep,[]);
        % contract UT*ST with M{itN}, to update M{itN}
        if itN > 1
            M{itN-1} = contract(M{itN-1},3,3,UT*diag(ST),2,1);
        end
        
        % update the Hamiltonian in effective basis
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
            % permute left<->right for Hs{itN} as well, to make use of updateLeft
            H2 = permute(Hs{itN},[3 2 1 4]); % right-bottom-left-top
            Hlr{itN+1} = updateLeft(Hlr{itN+2},3,T,H2,4,T);
        end
    end
    
    % display informaiton of the sweep
    str = ['Sweep #',sprintf('%i/%i',[2*itS-1,2*Nsweep]),...
        ' (right -> left) : Energy = ',sprintf('%.7g',Eiter(N,2*itS-1))];
    disptime(str);
    
    % left -> right
    for itN = (1:N)
        [Ceff,Eeff] = eigs_1site(Hlr{itN},Hs{itN},Hlr{itN+2},M{itN});
        
        Eiter(itN,2*itS) = Eeff;
        
        % update M{itN} and M{itN+1} by using Ceff, via SVD
        % decompose Ceff
        [M{itN},ST,VT] = svdTr(Ceff,3,[1 2],Nkeep,[]);
        % contract UT*ST with M{itN}, to update M{itN}
        if itN < N
            M{itN+1} = contract(diag(ST)*VT,2,2,M{itN+1},3,1);
        end
        
        % update the Hamiltonian in effective basis
        if itN == 1
            % "remove" the left leg (1st leg) of Hs{itN}, which is dummy,
            % by permuting to the last; MATLAB automatically suppresses the
            % trailing singleton dimensions
            H2 = permute(Hs{itN},[2 3 4 1]); % bottom-right-top (-left)
            Hlr{itN+1} = updateLeft([],[],M{itN},H2,3,M{itN});
        elseif itN == N
            Hlr{itN+1} = [];
        else
            Hlr{itN+1} = updateLeft(Hlr{itN},3,M{itN},Hs{itN},4,M{itN});
        end
    end
    
    % display informaiton of the sweep
    str = ['Sweep #',sprintf('%i/%i',[2*itS,2*Nsweep]),...
        ' (left -> right) : Energy = ',sprintf('%.7g',Eiter(N,2*itS))];
    disptime(str);

    % check energy convergence
    if itS > 1
        if abs(Eiter(N,2*itS) - Eiter(N,2*(itS-1))) < Econv
            break % if ((itS-1)th energy - (itS)th energy), stop DMRG sweep
        end
    end
end

E0 = Eiter(N,2*itS); % take the last value
    
toc2(tobj,'-v');

end