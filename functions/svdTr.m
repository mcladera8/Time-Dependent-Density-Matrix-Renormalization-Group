function [U,S,Vd,dw] = svdTr (T,rankT,idU,Nkeep,Stol,varargin)
% < Description >
%
% [U,S,Vd,dw] = svdTr (T,rankT,idU,Nkeep,Stol) % truncate by Nkeep and Stol
% [U,S,Vd,dw] = svdTr (T,rankT,idU,[],Stol) % truncate by Stol
% [U,S,Vd,dw] = svdTr (T,rankT,idU,Nkeep,[]) % truncate by Nkeep
% [U,S,Vd,dw] = svdTr (T,rankT,idU,[],[]) % not truncate (only the default tolerance Stol = 1e-8 is considered)
% [U,S,Vd,dw] = svdTr (.. ,'deg') % option to keep the degenerate states at trunctaion threshold
%
% Singular value decomposition of tensor such that T = U*diag(S)*Vd. (Note
% that it is not U*S*V' as in the MATLAB built-in function 'svd'.) If the
% truncation criterion is given, the tensors are truncated with respect to
% the largest singluar values.
%
% < Input >
% T : [tensor] Tensor.
% rankT : [number] Rank of T.
% idU : [integer vector] Indices of T to be associated with U. For example,
%       if rankT == 4 and idU == [1 3], the result U is rank-3 tensor whose
%       1st and 2nd legs correspond to the 1st and 3rd legs of T. The 3rd
%       leg of U is associated with the 1st leg of diag(S). And Vd is
%       rank-3 tensor whose 2nd and 3rd legs correspond to the 2nd and 4th
%       legs of T. Its 1st leg is associated with the 2nd leg of diag(S).
% Nkeep : [number] The number of singular values to keep.
%       (Default: Inf, i.e., no truncation)
% Stol : [number] Minimum magnitude of the singluar value to keep.
%       (Default: 1e-8, which is the square root of double precision 1e-16)
%
% < Option >
% 'deg' : If it is given, this funciton searches the largest gap in the log
%       of singular values, from the (Ntr)-th largetst singular value to
%       the ceil(Ntr*1.1)-th largest one. Use this option to avoid
%       unintended symmetry breaking by truncation.
%       (Default: off)
%
% < Output >
% U : [tensor] Tensor describing the left singular vectors. Its last leg
%       contracts with diag(S). The earlier legs are specified by input
%       idU; their order is determined by the ordering of idU.
% S : [vector] The column vector of singular values.
% Vd : [tensor] Tensor describing the right singular vectors. Its 1st leg
%       contracts with diag(S). The later legs conserve the order of the
%       legs of input T.
% dw : [tensor] Discarded weight (= sum of the square of the singular
%       values truncated).
%
% Written by S.Lee (May 22,2017)
% Updated by S.Lee (Jun.19,2017): Previosuly 'deg' option has been default;
%       changed to be optional.
% Updated by S.Lee (May 31,2019): Revised for Sose 2019.

% default truncation parameters
if isempty(Nkeep), Nkeep = Inf; end
if isempty(Stol), Stol = 1e-8; end
chkdeg = false;

% option
while ~isempty(varargin)
    switch varargin{1}
        case 'deg'
            chkdeg = true;
            varargin(1) = [];
        otherwise
            error('ERR: Unknown option.')
    end
end

Tdim = size(T); % dimensions of tensors

if rankT < numel(Tdim)
    error('ERR: Input ''rankT'' is smaller than the rank of other input ''T''.');
end
Tdim = [Tdim,ones(1,rankT-numel(Tdim))];

idTtot = (1:numel(Tdim));

% % % check integrity of inputs
if any(~any(bsxfun(@eq, idU(:), idTtot),2))
    error('ERR: Invalid index for tensor U (e.g. out of bound, non-integer).');
end
% % %

idV = idTtot; idV(idU) = [];
% reshape to matrix form
T2 = reshape(permute(T,[idU,idV]),[prod(Tdim(idU)) prod(Tdim(idV))]);
[U2,S2,V2] = svd(T2,'econ'); % SVD

S2 = diag(S2);

% number to be kept determined by Nkeep and by Stol
Ntr = [Nkeep,sum(S2>Stol)];
% choose stricter criterion; round up and compare with 1 just in case
Ntr = max(ceil(min(Ntr)),1); % keep at least one bond to maintain the tensor network structure

if chkdeg
    % indices of singluar values to find the largest gap
    ids_part = (max(Ntr,1):min(ceil(Ntr*1.1),numel(S2)));
    % in the above, max and min are given to be within range of indices
    if numel(ids_part) > 1
        [~,idmax] = min(diff(log(S2(ids_part)))); % find the largest gap (S is in decreasing order)
        Ntr = ids_part(idmax); % actual index of the largest kept energy eigenvalue
    end
end

dw = sum(S2(Ntr+1:end).^2);

S2 = S2(1:Ntr);
U2 = U2(:,(1:Ntr));
V2 = V2(:,(1:Ntr));

U = reshape(U2,[Tdim(idU) Ntr]);
S = S2;
Vd = reshape(V2',[Ntr Tdim(idV)]);

end