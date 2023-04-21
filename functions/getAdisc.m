function [odisc,Adisc] = getAdisc (Inrg,Op,Z,varargin)
% < Description >
%
% [odisc,Adisc] = getAdisc (Inrg,Op,Z [, option]); % for anti-commuting operators
% [odisc,Adisc] = getAdisc (Inrg,Op,[] [, option]); % for commuting operators
%
% Obtain the discrete spectral function, for given the operators at the
% impurity (site s00) and the NRG information. Using the Lehmann
% representation and the complete basis, we can obtain the spectral
% function as the collection of the delta functions with weights. The
% spectral weights are binned in the vector Adisc, where each bin of Adisc
% represents the narrow frequency interval centered at the corresponding
% element of odisc.
%
% < Input >
% Inrg : [struct] NRG information obtained after running NRG_1channel.
% Op : [tensor] Operator acting at the impurity (site s00).
% Z : [tensor] Fermionic sign operator. If the operator Op is commuting
%       (i.e., bosonic), put Z as empty [].
%
% < Option >
% 'emin', .. : [number] Minimum (in absolute value) frequency limit.
%               (Default: 1e-12)
% 'emax', .. : [number] Minimum (in absolute value) frequency limit.
%               (Default: 1e2)
% 'estep', .. : [number] Number of bins per decade (frequency value 1 ->
%               10).
%               (Default: 200)
%
% < Output >
% odisc : [vector] Logarithmic grid of frequency.
% Adisc : [vector] Discrete spectral weights corresponding to odisc.
%
% Written by S.Lee (May 22,2017)
% Updated by S.Lee (May 11,2019): Revised for SoSe 2019.
% Updated by S.Lee (Jun.20,2020): Revised for SoSe 2020.
% Updateb by J.Shim (Jul.21.2020): Revised for SoSe 2022.

tobj = tic2;

% default option values
emin = 1e-12;
emax = 1e2;
estep = 200;

% % parsing option
while ~isempty(varargin)
    switch varargin{1}
        case 'emin'
            emin = varargin{2};
            varargin(1:2) = [];
        case 'emax'
            emax = varargin{2};
            varargin(1:2) = [];
        case 'estep'
            estep = varargin{2};
            varargin(1:2) = [];
        otherwise
            error('ERR: Unknown option.');
    end
end
% % % 

if ~isempty(Z)
    opsign = 1; % anti-commutation
else
    opsign = -1; % commutation
end

logemin = log10(emin);
logemax = log10(emax);

odisc = 10.^(((logemin*estep):(logemax*estep))/estep);
Adisc = zeros(numel(odisc),2); % column 1: negative freq., col 2: positive freq.

if ~isempty(Z)
    disptime('Correlation function for anti-commuting op.');
else
    disptime('Correlation function for commuting op.');
end

% % Rerouting Z string: by doing this, we don't need to contract Z tensors
% when we update Op for longer chains. Refer to the appendix of A.
% Weichselbaum, Phys. Rev. B 86, 245124 (2012). (Especially Figs. 10 and
% 11)
if ~isempty(Z)
    Op = contract(Z,2,2,Op,3,1);
end

Oprev = updateLeft([],[],Inrg.AK{1},Op,3,Inrg.AK{1});
N = numel(Inrg.E0);

for itN = (2:N)
    Atmp = {Inrg.AD{itN},Inrg.AK{itN}};
    Rtmp = {diag(Inrg.RD{itN}),Inrg.RK{itN}};
    Etmp = {Inrg.ED{itN},Inrg.EK{itN}};
    
    for it1 = (1:2) % D(iscarded), K(ept)
        for it2 = (1:2) % D, K
            if ~isempty(Atmp{it1}) && ~isempty(Atmp{it2})
                % Onow: Operator in the current subspace
                Onow = updateLeft(Oprev,3,Atmp{it1},[],[],Atmp{it2});
                % ROp, OpR: contraction of density matrix and Onow
                ROp = contract(Rtmp{it1},2,2,Onow,3,1);
                OpR = contract(Onow,3,3,Rtmp{it2},2,1);

                if (it1 == 2) && (it2 == 2)
                    Oprev = Onow;
                else
                    % compute discrete spectral weights and add to the
                    % result 'Adisc'
                    Adisc = Adisc + getAdisc_1shell(Onow,ROp + opsign*OpR, ...
                        Etmp{it1},Etmp{it2},Inrg.EScale(itN), ...
                        numel(odisc),logemin,estep);
                end
            end
        end
    end
    
    disptime(['#',sprintf('%02i/%02i',[itN-1, N-1]),' : sum(Adisc) = ', ...
        sprintf('%.4g',sum(Adisc(:)))]); % show accumulative sum of the discrete weights
end

odisc = odisc(:);
odisc = [flipud(-odisc);odisc];
Adisc = [flipud(Adisc(:,1));Adisc(:,2)];

toc2(tobj,'-v');
chkmem;

end


function Adisc = getAdisc_1shell (Op1,Op2,E1,E2,EScale,nodisc,logemin,estep)
% compute spectral weights for given operators Op1 and Op2.

% permute the operator legs so that the operator-flavor leg is placed at
% the last
Op1 = permute(Op1,[1 3 2]);
Op2 = permute(Op2,[1 3 2]);

E21 = bsxfun(@minus,E2.',E1);
Op = sum(Op1.*conj(Op2),3);

E21 = E21(:);
Op = Op(:);

if numel(E21) ~= numel(Op)
    keyboard;
end

% indexing for frequency
ids = round((log10(abs(E21))+log10(EScale)-logemin)*estep)+1;
ids(ids<1) = 1;
ids(ids>nodisc) = nodisc;

% indexing for sign
sgn = (E21 >= 0)+1;

Adisc = accumarray([ids,sgn],Op,[nodisc,2]);

end