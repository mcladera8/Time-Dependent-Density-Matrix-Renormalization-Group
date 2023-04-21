function [M,S,dw] = canonForm (M,id,varargin)
% < Description >
%
% [M,S,dw] = canonForm (M,id [,Nkeep] [,'-k']);
%
% Obtain the canonical forms of MPS. It brings the tensors M{1}, ..., M{id}
% into the left-canonical form and the others M{id+1}, ..., M{end} into the
% right-canonical form.
%
% < Input >
% M : [cell array] MPS of length numel(M). Each cell element is a rank-3
%       tensor, where the first, second, and third dimensions are
%       associated with left, bottom (i.e., local), and right legs,
%       respectively.
% id : [integer] Index for the bond connecting the tensors M{id} and
%       M{id+1}. With respect to the bond, the tensors to the left
%       (right) are brought into the left-(right-)canonical form. If id ==
%       0, the whole MPS will be in the right-canonical form.
%
% < Option >
% Nkeep : [integer] Maximum bond dimension. That is, only Nkeep the
%       singular values and their associated singular vectors are kept at
%       each iteration.
%       (Default: Inf)
% '-k' : If set, the very small singular values (smaller than eps) and
%       their corresponding singular vectors are not truncated.
%       (Default: not set, i.e., truncate very small singular values)
%
% < Output >
% M : [cell array] Left-, right-, or bond-canonical form from input M,
%       depending on id, as follows:
%       * id == 0: right-canonical form
%       * id == numel(M): left-canonical form
%       * otherwise: bond-canonical form
% S : [column vector] Singular values at the bond between M{id} and M{id+1}. 
% dw : [column vector] Vector of length numel(M)-1. dw(n) means the
%       discarded weight (i.e., the sum of the square of the singular  
%       values that are discarded) at the bond between M{n} and M{n+1}.
%
% Written by S.Lee (Apr.30,2019)
% Revised by S.Lee (Jul.21,2021): Added an option '-k'.

% % default option values
Nkeep = Inf; % keep all
is2keep0 = false; % false: truncate very small singular values and their corresponding singular vectors

% % parsing option
while ~isempty(varargin)
    if isequal(varargin{1},'-k')
        is2keep0 = true;
    else   
        Nkeep = varargin{1}(1);
    end
    varargin(1) = [];
end
% % % %

% % check the integrity of input
if (numel(id) ~= 1) || (round(id) ~= id)
    error('ERR: 2nd input ''id'' needs to be a single integer.');
elseif (id < 0) || (id > numel(M))
    error('ERR: the 2nd input ''id'' needs to be in a range (0:numel(M))');
elseif size(M{1},1) ~= 1
    error('ERR: the first dimension (= left leg) of M{1} should be of size 1.');
elseif size(M{end},3) ~= 1
    error('ERR: the third dimension (= right leg) of M{end} should be of size 1.');
elseif (Nkeep <= 1) || isnan(Nkeep)
    error('ERR: Option ''Nkeep'' should be positive integer larger than 1.');
end
% % % %

dw = zeros(numel(M)-1,1); % discarded weights

% % Bring the left part of MPS into the left-canonical form
for it = (1:id)
    % % % TODO (start) % % %
    
    % reshape M{it} and SVD
    T = M{it};
    T = reshape(T,[size(T,1)*size(T,2),size(T,3)]);
    [U,S,V] = svd(T,'econ');
    
    % remove the components of U, S, V associated with singular values
    % smaller than the double precision noise
    Svec = diag(S); % vector of singular values
    if ~is2keep0
        ok = (Svec < eps);
        U(:,ok) = [];
        Svec(ok) = [];
        V(:,ok) = [];
    end
    
    % truncate singular values/vectors; keep up to Nkeep. Truncation at the
    % bond between M{id} and M{id+1} is performed later.
    if ~isinf(Nkeep) && (it < id)
        nk = min([numel(Svec);Nkeep]); % actual number of singular values/vectors to keep
        dw(it) = dw(it) + sum(Svec(nk+1:end).^2); % discarded weights
        U = U(:,(1:nk));
        V = V(:,(1:nk));
        Svec = Svec(1:nk);
    end
    
    S = diag(Svec); % return to square matrix
    
    % reshape U into rank-3 tensor, and replace M{it} with it
    M{it} = reshape(U,[size(U,1)/size(M{it},2),size(M{it},2),size(U,2)]);
    
    if it < id
        % contract S and V' with M{it+1}
        M{it+1} = contract(S*V',2,2,M{it+1},3,1);
    else
        % R1: tensor which is the leftover after transforming the left
        %   part. It will be contracted with the counterpart R2 which is
        %   the leftover after transforming the right part. Then R1*R2 will
        %   be SVD-ed and its left/right singular vectors will be
        %   contracted with the neighbouring M-tensors.
        R1 = S*V';
    end
    
    % % % TODO (end) % % %
end

% % In case of fully right-canonical form; the above for-loop is not executed
if id == 0
    R1 = 1;
end
    
% % Bring the right part into the right-canonical form
for it = (numel(M):-1:id+1)
    % % % TODO (start) % % %
    
    % reshape M{it} and SVD
    T = M{it};
    T = reshape(T,[size(T,1),size(T,2)*size(T,3)]);
    [U,S,V] = svd(T,'econ');
    
    % remove the components of U, S, V associated with singular values
    % smaller than the double precision noise
    Svec = diag(S); % vector of singular values
    if ~is2keep0
        ok = (Svec < eps);
        U(:,ok) = [];
        Svec(ok) = [];
        V(:,ok) = [];
    end
    
    % truncate singular values/vectors; keep up to Nkeep. Truncation at the
    % bond between M{id} and M{id+1} is performed later.
    if ~isinf(Nkeep) && (it > (id+1))
        nk = min([numel(Svec);Nkeep]); % actual number of singular values/vectors to keep
        dw(it-1) = dw(it-1) + sum(Svec(nk+1:end).^2); % discarded weights
        U = U(:,(1:nk));
        V = V(:,(1:nk));
        Svec = Svec(1:nk);
    end
    
    S = diag(Svec); % return to square matrix
    
    % reshape V' into rank-3 tensor, replace M{it} with it
    M{it} = reshape(V',[size(V,2),size(M{it},2),size(V,1)/size(M{it},2)]);
    
    if it > (id+1)
        % contract U and S with M{it-1}
        M{it-1} = contract(M{it-1},3,3,U*S,2,1);
    else
        % R2: tensor which is the leftover after transforming the right
        %   part. See the description of R1 above.
        R2 = U*S;
    end
    
    % % % TODO (end) % % %
end

% % In case of fully left-canonical form; the above for-loop is not executed
if id == numel(M)
    R2 = 1;
end

% % SVD of R1*R2, and contract the left/right singular vectors to the tensors
[U,S,V] = svd(R1*R2,'econ');
S = diag(S); % vector of singular values

% remove the components of U, S, V associated with singular values
% smaller than the double precision noise
if ~is2keep0
    ok = (S < eps);
    U(:,ok) = [];
    S(ok) = [];
    V(:,ok) = [];
end

% truncate singular values/vectors; keep up to Nkeep. At the leftmost and
% rightmost legs (dummy legs), there should be no truncation, since they
% are already of size 1.
if ~isinf(Nkeep) && (id > 0) && (id < numel(M))
    % % % TODO (start) % % %
    
    nk = min([numel(S);Nkeep]); % actual number of singular values/vectors to keep
    dw(id) = dw(id) + sum(S(nk+1:end).^2); % discarded weights
    U = U(:,(1:nk));
    V = V(:,(1:nk));
    S = S(1:nk);
    
    % % % TODO (end) % % %
end

if id == 0 % fully right-canonical form
    % U is a single number which serves as the overall phase factor to the
    % total many-site state. So we can pass over U to V'.
    M{1} = contract(U*V',2,2,M{1},3,1);
elseif id == numel(M) % fully left-canonical form
    % V' is a single number which serves as the overall phase factor to the
    % total many-site state. So we can pass over V' to U.
    M{end} = contract(M{end},3,3,U*V',2,1);
else
    M{id} = contract(M{id},3,3,U,2,1);
    M{id+1} = contract(V',2,2,M{id+1},3,1);
end

end