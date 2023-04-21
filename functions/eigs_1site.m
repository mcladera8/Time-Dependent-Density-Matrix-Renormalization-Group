function [Ceff,Eeff] = eigs_1site(Hleft,Hloc,Hright,Cinit,varargin)
% < Description >
%
% [Ceff,Eeff] = eigs_1site(Hleft,Hloc,Hright,Cinit [, option])
%
% Obtain the ground state and its energy for the effective 
% Hamiltonian for the site-canonical MPS, by using the Lanczos 
% method.
%
% < Input >
% Hleft, Hloc, Hright: [tensors] The Hamiltonian for the left, site
%	 and right parts of the chain. They form the effective 
%	 Hamiltonian in the site-canonical basis.
% Cinit : [tensor] Ket tensor at a lattice site. It becomes an 
%	initial vector for the Lanczos method.
%
% The input tensors can be visualized as follows:
% (numbers are the order of legs, * means contraction)
%
%
%      	    1 -->-[ Cinit ]-<-- 3
%                     |
%                     ^ 2
%                     |
%
%
%     /--->- 3        | 4        3 -<---\
%     |               ^                 |
%     |    2     1    |    3     2      |
%   Hleft-->- * -->- Hloc->-- * ->-- Hright
%     |               |                 |
%     |               ^                 |
%     \---<- 1        | 2        1 ->---/
%
% < Option >
% 'N', .. : [numeric] Maximum number of Lanczos vectors (in 
%	addition to those given by Cinit) to be considered 
%	for the Krylov subspace.
%       (Default: 5)
% 'minH', .. : [numeric] Minimum absolute value of the 1st diagonal 
%	(i.e., superdiagonal) element of the Hamiltonian in the 
%	Krylov subspace. If a 1st-diagonal element whose absolute 
%	value is smaller than minH is encountered, the iteration 
%	stops. Then the ground-state vector and energy is obtained 
%	from the tridiagonal matrix constructed so far.
%       (Default: 1e-10)
%
% < Output >
% Ceff : [tensor] A ket tensor as the ground state of the 
%	effective Hamiltonian.
% Eeff : [numeric] The energy eigenvalue corresponding to Ceff.
% Written by S.Lee (May 31,2017)
% Documentation updated by S.Lee (Jun.8,2017)
% Updated by S.Lee (May 28,2019): Revised for SoSe 2019.
% Updated by S.Lee (May 23,2020): Revised for SoSe 2020.
% Updated by S.Lee (Jun.08,2020): Typo fixed.
% Updated by J.Shim (May 25,2022): Revised for SoSe 2022.

% default parameter
N = 5;
minH = 1e-10;

% parsing option
while ~isempty(varargin)
    switch varargin{1}
        case 'N'
            N = varargin{2};
            varargin(1:2) = [];
        case 'minH'
            minH = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: check input!');
    end
end

% size of ket tensor
Csz = [size(Cinit,1),size(Cinit,2),size(Cinit,3)];

% initialize Cinit
Cinit = Cinit/norm(Cinit(:)); % normalize Cinit

% Krylov vectors (vectorized tensors)
Ckr = zeros(numel(Cinit),N+1);
Ckr(:,1) = Cinit(:);

% In the Krylov basis, the Hamiltonian becomes tridiagonal
ff = zeros(N,1); % 1st diagonal
gg = zeros(N+1,1); % main diagonal

for itN = (1:(N+1))
    % contract Hamiltonian with ket tensor
    Ctmp = eigs_1site_HC(Hleft,Hloc,Hright,reshape(Ckr(:,itN),Csz));
    Ctmp = Ctmp(:); % vectorize
    
    gg(itN) = Ckr(:,itN)'*Ctmp; % diagonal element; "on-site energy"
    
    if itN < (N+1)
        % orthogonalize Atmp w.r.t. the previous ket tensors
        Ctmp = Ctmp - Ckr(:,(1:itN))*(Ckr(:,(1:itN))'*Ctmp);
        % twice, to reduce numerical noise
        Ctmp = Ctmp - Ckr(:,(1:itN))*(Ckr(:,(1:itN))'*Ctmp);
        
        % norm
        ff(itN) = norm(Ctmp);
        
        if ff(itN) > minH
            Ckr(:,itN+1) = Ctmp/ff(itN);
        else
            % stop iteration; truncate ff, gg
            ff(itN:end) = [];
            gg(itN+1:end) = [];
            Ckr(:,(itN+1):end) = [];
            break;
        end
    end
end

% Hamiltonian in the Krylov basis
Hkr = diag(ff,1);
Hkr = Hkr + Hkr' + diag(gg);
[Vkr,Ekr] = eig((Hkr+Hkr')/2);
[~,minid] = min(diag(Ekr));

% ground state
Ceff = Ckr*Vkr(:,minid);
Ceff = Ceff/norm(Ceff); % normalize
Ceff = reshape(Ceff,Csz); % reshape to rank-3 tensor

% ground-state energy; measure again
Ctmp = eigs_1site_HC(Hleft,Hloc,Hright,Ceff);
Eeff = Ceff(:)'*Ctmp(:);

end

function Cout = eigs_1site_HC(Hleft,Hloc,Hright,Cin)
% < Description >
%
% Cout = eigs_1site_HC(Hleft,Hloc,Hright,Cin)
%
% Apply the effective Hamitonian for the site-canonical MPS.
% 
% < Input >
% Hleft, Hloc, Hright: [tensors] Refer to the description of the variables
%       with the same names, in 'DMRG_1site_eigs'.
% Cin : [tensor] A ket tensor at a lattice site, to be applied by the
%       effective Hamiltonian.
%
% < Output >
% Cout : [tensor] A ket tensor at a lattice site, after the application of
%   the effective Hamiltonian to Cin.
%
% Written by S.Lee (May 23,2020)
% Updated by S.Lee (May 27,2020): Minor change.
% Updated by J.Shim (May 25,2022): Revised for SoSe 2022.

% set empty tensors as 1, for convenience
if isempty(Hleft)
    Hleft = 1;
end
if isempty(Hright)
    Hright = 1;
end

Cout = contract(Hleft,3,3,Cin,3,1);
Cout = contract(Cout,4,[2 3],Hloc,4,[1 4]);
Cout = contract(Cout,4,[2 4],Hright,3,[3 2]);
    
end