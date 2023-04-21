function Inrg = getRhoFDM (Inrg,T)
% < Description >
%
% Inrg = getRhoFDM (Inrg,T)
%
% Construct the full density matrix (FDM) in the basis of both discarded
% and kept states, for given temperature T.
%
% < Input >
% Inrg : [struct] NRG information obtained after running NRG_1channel.
% T : [number] Temperature. Here we set \hbar = k_B = 1.
%
% < Ouput >
% Inrg : [struct] NRG result. It keeps the result of NRG_1channel. In
%       addition to the result, this function adds two more fields to Inrg:
%   .RD, .RK : [cell] Full density matrix in the discarded and kept state
%       basis, respectively. Each cell element Inrg.RD{n} is a column
%       vector whose elements are the density matrix elements associated
%       with the discarded energy eigenstates at the iteration n-1. (Note
%       that s00 for n = 1 is for the iteration diagonalizing K00 basis.)
%       Inrg.RK{n} is a matrix in the basis of the kept energy eigenstates
%       at the iteration n-1.
%
% Written by S.Lee (May 22,2017)
% Updated by S.Lee (May 12,2019): Revised for SoSe 2019.
% Updated by S.Lee (Jun.20,2020): Revised for SoSe 2020.
% Updated by J.Shim (Jul.07.2022): Revised for SoSe 2022.

tobj = tic2;
disptime(['Construct full density matrix @ T = ',sprintf('%.4g',T),' ...']);

L = numel(Inrg.E0);

% extract the local space dimension from ket tensors
locdim = zeros(L,1);
for itL = (1:L)
    if ~isempty(Inrg.AK{itL})
        locdim(itL) = size(Inrg.AK{itL},2);
    else
        locdim(itL) = size(Inrg.AD{itL},2);
    end
end

% the shift of energy in each shell measured from the lowest-energy of the
% last iteration
E0r = [Inrg.EScale(2:end).*Inrg.E0(2:end),0];
E0r = fliplr(cumsum(fliplr(E0r)));

RD = cell(1,L); % FDM in the discarded state basis; row vector
RK = cell(1,L); % FDM in the kept state basis; matrix

RDsum = zeros(1,L); % sum of Boltzmann weights

% obtain the Boltzamann weights
for itL = (1:L)
    % Obtain the column vector RD{itL} whose elements are the Boltzmann
    % weights
    RD{itL} = exp(-(Inrg.ED{itL}*Inrg.EScale(itL)-E0r(itL))/T)*prod(locdim(itL+1:end));
    RDsum(itL) = sum(RD{itL});
end

RDsum = sum(RDsum);

% normalize the Boltzmann weights to get the elements of the density matrix
% in the discarded basis
for itL = (1:L)
    RD{itL} = RD{itL}/RDsum;
end

% update the FDM in the kept basis
for itL = (L:-1:2)
    % Construct RK{itL-1} as the sum of RD{itL} and RK{itL}, with the local
    % Hilbert space for the site s(itL-1). (Note that s00 for itL = 1, s01
    % for itL = 2, etc.)
    % Hint: one may utilize updateLeft, after permuting the legs of
    % Inrg.AD{itL} and Inrg.AK{itL}.
    
    % to utilize updateLeft for contracting AD with RD, and AK with RK,
    % permute the legs of AD and AK (left <-> right)
    AD2 = permute(Inrg.AD{itL},[3 2 1]);
    AK2 = permute(Inrg.AK{itL},[3 2 1]);
    RK{itL-1} = updateLeft(diag(RD{itL}),2,AD2,[],[],AD2) ...
        + updateLeft(RK{itL},2,AK2,[],[],AK2);
    % NOTE: AK and AD are in left-canonical form, not right-canonical.
end

if sum(RD{end}) > 1e-2
    disptime('WRN: sum(Inrg.RD{end}) > 1e-2 ; chain length is not enough');
end

Inrg.RK = RK;
Inrg.RD = RD;

toc2(tobj,'-v');

end