%% All COMMENTS go here

%-------- Intro and doLinDisc

%{
% N equidistant intervals Ik of size ∆e. 
% e_k = min{Ik} + ∆e/2
clear
% Create the semicircular bath
N = 40;
ozin = linspace(-1,1,N+1);
deltaE = ozin(2)-ozin(1);
for i=(1:numel(ozin))
    RhoV2in(i)=(1/(2*pi()))*sqrt(1-(ozin(i))^2);
end
%figure
%clc
%plot(ozin,RhoV2in)
for i=(1:N)
    estar(i) = min(ozin(i),ozin(i+1)) + deltaE/2;
end
%figure
%plot(estar,'+')

for i=(1:N)
    a = [ozin(i) ozin(i+1)];
    b = [RhoV2in(i) RhoV2in(i+1)];
    vstar(i) = sqrt(trapz(a,b));
end
%figure
%plot((1:N-1),vstar,'ro')

estarimp = [fliplr(estar(N/2+1:N)) 0 fliplr(estar(1:N/2))];
vstarimp = [vstar(1:N/2) 0 vstar(N/2+1:N)];
Hstar = diag(estarimp);
Hstar(:,N/2+1) = vstarimp;
Hstar(N/2+1,:) = vstarimp;


U = zeros(size(Hstar,1),1);
U(N/2 +1,1) = 1;

% Lanczos tridiagonalization
ff = zeros(N,1); % hopping amplutudes; corresponds to the super diagonal 
                  % (diag(..,+1) or diag(..,-1)) in the tridiagonal matrix
gg = zeros(N,1); % hopping amplutudes; corresponds to the main diagonal 
                  % (diag(..)) in the tridiagonal matrix


for itN = (1:N)
    v = Hstar*U(:,itN);
    v = v-U*(U'*v);
    v = v-U*(U'*v); % twice for numerical reason
    ff(itN) = norm(v);

    if (itN < N) && (ff(itN) > 0)
        U(:,itN+1) = v/ff(itN);
        gg(itN) = U(:,itN+1)'*Hstar*U(:,itN+1);
    elseif (itN < N) && (ff(itN) <= 0)
        ff = ff(1:itN);
        gg = gg(1:itN);
        break
    end
end


figure
set(gcf,'units','points');
set(gcf,'position',[0 0 450 600]);
set(gcf,'paperunits', 'points');
set(gcf,'papersize', [450 600]);
set(gcf,'paperpositionmode', 'manual');
set(gcf,'paperposition', [0 0 450 600]);
set(gca,'units','points');
set(gca,'innerposition',[50,50,350,500]);
plot((1:N),ff,'+','linewidth',2,'MarkerSize',10);
hold on
plot((1:N),6*vstar,'x','linewidth',2,'MarkerSize',10);
hold off
xlabel('$\ell$','Interpreter','latex','fontsize',18);
xlim([0 40.5])
set(gca,'FontSize',18);
set(gca,'yminortick','on');
set(gca,'ticklength',[0.04 0.1]);
legend('$t_{chain}$',...
    '$6\cdot t_{star}$', ...
    'Interpreter','latex','fontsize',18,'location','south');

figure
set(gcf,'units','points');
set(gcf,'position',[0 0 450 600]);
set(gcf,'paperunits', 'points');
set(gcf,'papersize', [450 600]);
set(gcf,'paperpositionmode', 'manual');
set(gcf,'paperposition', [0 0 450 600]);
set(gca,'units','points');
set(gca,'innerposition',[50,50,350,500]);
plot((1:N),gg(1:N),'+','linewidth',2,'MarkerSize',10);
hold on
plot((1:N),estar,'x','linewidth',2,'MarkerSize',10);
hold off
xlabel('$\ell$','Interpreter','latex','fontsize',18);
xlim([0 40.5])
set(gca,'FontSize',18);
set(gca,'yminortick','on');
set(gca,'ticklength',[0.04 0.1]);
legend('$\varepsilon_{chain}$',...
    '$\varepsilon_{star}$', ...
    'Interpreter','latex','fontsize',18,'location','south');



%% For the tDMRG function in case it was as I did it before 
This is definetly not ok
    Gr = 0;
    if mod(it1,3) == 0
        for i=(1:N)
            Gr = updateLeft([],[],Mf{i},[],[],conj(Mr{i}))*exp(1i*Econv*2*dt*it1/3);
        end
        Green(round(it1/3)) = Gr;
    end
%}

%-------- Eigendescomp Chain and Star

%{

[~,GreenChain] = eig((Hchain+Hchain')/2);
[~,GreenStar] = eig((Hstar+Hstar')/2);
GreenChain = diag(GreenChain);
GreenStar = diag(GreenStar);
%}

%-------- MPO and Iterative Diag

%{
% MPO for Hff: START %
% As stated in the paper by Daniel Bauernfeind et.al., we'll be using a
% local Hilbert of dimension 2 (empty |0> or occupied |1>)
 W = cell(N+1,1); % MPO W
 % First site
 W{1} = zeros(1,2,4,2); % ordering: left bottom right top
 W{1}(1,:,2,:) = ff(1)*squeeze(F)'; % t_l*creation
 W{1}(1,:,3,:) = ff(1)'*squeeze(F); % t_l*annihilation
 W{1}(1,:,4,:) = I;
    % Last site
 W{N+1} = zeros(4,2,1,2); % ordering: left bottom right top
 W{N+1}(1,:,1,:) = I;
 W{N+1}(2,:,1,:) = squeeze(F); % annihilation
 W{N+1}(3,:,1,:) = squeeze(F)'; % creation
 % Other sites
 for itL2 = (2:N)
        W{itL2} = zeros(4,2,4,2); % ordering: left bottom right top
        W{itL2}(1,:,1,:) = I;
        W{itL2}(2,:,1,:) = squeeze(F); % annihilation
        W{itL2}(3,:,1,:) = squeeze(F)'; % creation
        W{itL2}(4,:,2,:) = ff(itL2)*squeeze(F)'; % t_l*creation
        W{itL2}(4,:,3,:) = ff(itL2)'*squeeze(F); % t_l*annihilation
        W{itL2}(4,:,4,:) = I;
 end
%}
%{
% GS MPS for Hff: START %
% In order to obtain the first gsMPS we do iterative diagonalization so
% that then variationally (via DMRG) we obtain the final Ground State.
tol = Nkeep*100*eps; % numerical tolerance for degeneracy
H0 = I*0; % 1st site Hamiltonian
A0 = getIdentity(1,2,I,2); % 1st leg: dummy
MPS = cell(N+1,1);
Hnow = H0;
[V,D] = eig((Hnow+Hnow')/2);
MPS{1} = contract(A0,3,3,V,2,1);
Hprev = D;
for itL2 = (2:N+1)
        % Fermion aniihilation operator at the current site
        Fprev = updateLeft([],[],MPS{itL2-1},F,3,MPS{itL2-1});
        Anow = getIdentity(Hprev,2,I,2);
        Hnow = updateLeft(Hprev,2,Anow,[],[],Anow);
        % hopping terms
        Hhop = ff(itL2-1)*updateLeft(Fprev,3,Anow,permute(F,[3 2 1]),3, ...
            Anow);
        Hhop = Hhop + Hhop';
        Hnow = Hnow + Hhop;
        [V,D] = eig((Hnow+Hnow')/2);
        [D,ids] = sort(diag(D),'ascend');
        V = V(:,ids);
        % truncation threshold for energy
        Etr = D(min([numel(D);Nkeep]));
        oks = (D <= (Etr + tol));
        if itL2 < N+1
            MPS{itL2} = contract(Anow,3,3,V(:,oks),2,1);
        elseif itL2 == N+1
            % Choose GS at the last step
            MPS{itL2} = contract(Anow,3,3,V(:,1),2,1);
        end
        Hprev = diag(D(oks));
%         disptime(['#',sprintf('%02i/%02i',[itL2,L]),' : ', ...
%             'NK=',sprintf('%i/%i',[size(MPS{itL2},3),size(Hnow,2)])]);
end
% GS MPS for Hff: END %
%}

%-------- MPS

%{
M{1} = contract(M{1},3,2,Fdag,2,2,[1 3 2]);
for i=(2:(numel(M)))
    M{i} = contract(M{i},3,2,Z,2,2,[1 3 2]);
end
%}
%{    
% Creation operator MPO for only site 1 
    CreI = cell(Nsite,1);
    % First Site
    Crei{1} = zeros(1,2,2,2);
    Crei{1}(1,:,1,:) = Fdag;
    Crei{1}(1,:,2,:) = eye(2);

    % Last Site
    Crei{Nsite} = zeros(2,2,1,2);
    Crei{Nsite}(1,:,1,:) = Z;

    % Intermidiate Sites
    for iN = (2:(Nsite-1))
        Crei{iN} = zeros(2,2,2,2);
        Crei{iN}(1,:,1,:) = Z;
        Crei{iN}(2,:,2,:) = eye(2);
    end
%}

%-------- Entropy figures

%{
figure
plot((1:(Nsite-1)),EE(75,:))


figure
set(gcf,'units','points');
set(gcf,'position',[0 0 375 310]);
set(gcf,'paperunits', 'points');
set(gcf,'papersize', [975 310]);
set(gcf,'paperpositionmode', 'manual');
set(gcf,'paperposition', [0 0 375 310]);
set(gca,'units','points');
imagesc((1:N-1),Tstep,EE);
xlabel('site $\ell$','Interpreter','latex','fontsize',18);
ylabel('time $t$','Interpreter','latex','fontsize',18);
colormap(turbo)
colorbar('eastoutside');
set(gca,'innerposition',[70,50,220,220]);
set(gca,'FontSize',18);
title('Entanglement Entropy','interpreter','latex');



%}

%-------- Green Exact

%{
[~,GreenChain]=eig((Hchain+Hchain')/2);
GreenChain = diag(GreenChain);

GreenChainEVO = zeros(length(Tstep),1);
for i=(1:length(Tstep))
    GreenChainEVO(i) = GreenChain(1)*exp(1i*E0*2*dt*i);
end
figure
plot(2*Tstep(1:20),real(GreenChainEVO(1:20)))

% -------
[Msr,EE2,dw,~] = tDMRG_star(M,Hs,Nkeep,-dt,-tmax,E0);
%}
%{
[D,GreenChain]=eig((Hchain+Hchain')/2);
GreenChain = diag(GreenChain);

GreenChainEVO = zeros(length(Tstep),1);
for i=(1:length(Tstep))
    GreenChainEVO(i) = GreenChain(1)*exp(1i*E0*2*dt*i);
end

figure
plot(2*Tstep,real(GreenChainEVO))
%}
%{
dt = 0.05;
xt = 0:dt:18;
[U1,E1] = eig((Hchain+Hchain')/2);
Ud = U1';
Ed = diag(E1);
Cn2 = zeros(Nsite,1);  % |C_n|^2
Crx = cell(Nsite,1);  % Creation MPOs.
Anx = cell(Nsite,1);  % Annihilation MPOs.
sm = squeeze(F);
sp = squeeze(F)';
sz = squeeze(Z);
for itk = 1:(Nsite)
    % First: build annihilation operator MPO.
    Anx{1} = zeros(1,2,2,2);
    Anx{1}(1,:,1,:) = Ud(itk,1)*sm;
    Anx{1}(1,:,2,:) = eye(2);
    for itj = 2:(Nsite-1)
        Anx{itj} = zeros(2,2,2,2);
        Anx{itj}(1,:,1,:) = sz;
        Anx{itj}(2,:,1,:) = Ud(itk,itj)*sm;
        Anx{itj}(2,:,2,:) = eye(2);
    end
    Anx{2*N+2} = zeros(2,2,1,2);
    Anx{2*N+2}(1,:,1,:) = sz;
    Anx{2*N+2}(2,:,1,:) = Ud(itk,2*N+2)*sm;
    % Second: build creation operator MPO.
    Crx{1} = zeros(1,2,2,2);
    Crx{1}(1,:,1,:) = U1(1,itk)*sp;
    Crx{1}(1,:,2,:) = eye(2);
    for iti = 2:(Nsite-1)
        Crx{iti} = zeros(2,2,2,2);
        Crx{iti}(1,:,1,:) = sz;
        Crx{iti}(2,:,1,:) = U1(iti,itk)*sp;
        Crx{iti}(2,:,2,:) = eye(2);
    end
    Crx{Nsite} = zeros(2,2,1,2);
    Crx{Nsite}(1,:,1,:) = sz;
    Crx{Nsite}(2,:,1,:) = U1(2*N+2,itk)*sp;
    % Apply the number operator to Psi0.
    % % contract the left half.
    cl = contract(Anx{1},4,4,M{1},3,2);
    cl = contract(Crx{1},4,4,cl,5,2);
    cl = contract(conj(M{1}),3,2,cl,7,2);
    cl = permute(cl,[2,4,6,8,1,3,5,7]);
    for j = 2:(N+1)
        cl = contract(cl,4,4,M{j},3,1);
        cl = contract(Anx{j},4,[1,4],cl,5,[3,4]);
        cl = contract(Crx{j},4,[1,4],cl,5,[4,1]);
        cl = contract(conj(M{j}),3,[1,2],cl,5,[4,1]);
    end
    if numel(size(cl)) < 3
        cl = reshape(cl,[1,size(a)]);
    end
    % % contract the right half.
    cr = contract(Anx{2*N+2},4,4,M{2*N+2},2,2);
    cr = contract(Crx{2*N+2},4,4,cr,4,2);
    cr = contract(conj(M{2*N+2}),2,2,cr,6,2);
    cr = permute(cr,[1,2,4,6,3,5]);
    for j = (2*N+1):-1:(N+2)
        cr = contract(cr,4,4,M{j},3,3);
        cr = contract(Anx{j},4,[3,4],cr,5,[3,5]);
        cr = contract(Crx{j},4,[3,4],cr,5,[4,2]);
        cr = contract(conj(M{j}),3,[2,3],cr,5,[2,4]);
    end
    % % combine left and right halves.
    Cn2(itk) = contract(cl,4,[1,2,3,4],cr,4,[1,2,3,4]);
end

for itT = 1:numel(xt)
    Gexact(itT) = sum(Cn2.*exp(-1i*(Ed-E0)*xt(itT)));
end
figure
plot(xt,real(Gexact))
%}
%{
M = canonForm(M,Nsite,'-k');


dt = 0.1;
xt = 0:dt:20;
[U1,E1] = eig((Hchain+Hchain')/2);
Ud = U1';
Ed = diag(E1);
Cn2 = zeros(Nsite,1);  % |C_n|^2
Crx = cell(Nsite,1);  % Creation MPOs.
Anx = cell(Nsite,1);  % Annihilation MPOs.
sm = squeeze(F);
sp = squeeze(F)';
sz = squeeze(Z);
for itk = 1:(Nsite)
    % First: build annihilation operator MPO.
    Anx{1} = zeros(1,2,2,2);
    Anx{1}(1,:,1,:) = Ud(itk,1)*sm;
    Anx{1}(1,:,2,:) = eye(2);
    for itj = 2:(Nsite-1)
        Anx{itj} = zeros(2,2,2,2);
        Anx{itj}(1,:,1,:) = sz;
        Anx{itj}(2,:,1,:) = Ud(itk,itj)*sm;
        Anx{itj}(2,:,2,:) = eye(2);
    end
    Anx{Nsite} = zeros(2,2,1,2);
    Anx{Nsite}(1,:,1,:) = sz;
    Anx{Nsite}(2,:,1,:) = Ud(itk,Nsite)*sm;
    % Second: build creation operator MPO.
    Crx{1} = zeros(1,2,2,2);
    Crx{1}(1,:,1,:) = U1(1,itk)*sp;
    Crx{1}(1,:,2,:) = eye(2);
    for iti = 2:(Nsite-1)
        Crx{iti} = zeros(2,2,2,2);
        Crx{iti}(1,:,1,:) = sz;
        Crx{iti}(2,:,1,:) = U1(iti,itk)*sp;
        Crx{iti}(2,:,2,:) = eye(2);
    end
    Crx{Nsite} = zeros(2,2,1,2);
    Crx{Nsite}(1,:,1,:) = sz;
    Crx{Nsite}(2,:,1,:) = U1(Nsite,itk)*sp;
    % Apply the number operator to Psi0.
        % % contract the left half.
    cl = contract(Anx{1},4,4,M{1},3,2);
    cl = contract(Crx{1},4,4,cl,5,2);
    cl = contract(conj(M{1}),3,2,cl,7,2);
    cl = permute(cl,[2,4,6,8,1,3,5,7]);
    for j = 2:(Nsite)
        cl = contract(cl,4,4,M{j},3,1);
        cl = contract(Anx{j},4,[1,4],cl,5,[3,4]);
        cl = contract(Crx{j},4,[1,4],cl,5,[4,1]);
        cl = contract(conj(M{j}),3,[1,2],cl,5,[4,1]);
    end
    Cn2(itk) = cl;
    % % combine left and right halves.
    %Cn2(itk) = contract(cl,4,[1,2,3,4],cr,4,[1,2,3,4]);
end
sum(Cn2)
for itT = 1:numel(xt)
    Gexact(itT) = sum(Cn2.*exp(-1i*(Ed-E0)*xt(itT)));
end
figure
plot(xt,real(Gexact))
%}

%-------- tDMRG function

%{
function Ovals = tDMRG_expVal(M,O,isleft)
%{
 Calculates the greater Green's function.
 < Input >
 Mf(Mr) : [cell] Input MPS's.
 Egs : [numerical] 
 isleft : [logical] If true, it means that the MPS M is in left-canonical
       form. Otherwise, right-canonical form.
 dt : [numerical] Time-Step
 < Output >
 Ovals : [vector] Green function at the desired time step.
%}

N = length(M);

Ovals = zeros(N,1);

MM = []; % contraction of bra/ket tensors
if isleft ~= 1 	% right-normalized
    for itN = (N:-1:1)
        T = permute(M{itN},[3 2 1]); % permute left<->right to use uLeft                                        
        Ovals(itN) = trace(updateLeft(MM,2,T,O,2,T));
        MM = updateLeft(MM,2,T,[],[],T);
    end
else 	% left-normalized
    for itN = (1:N)
        Ovals(itN) = trace(updateLeft(MM,2,M{itN},O,2,M{itN}));
        MM = updateLeft(MM,2,M{itN},[],[],M{itN});
    end
end


end
%}
%{
%*exp(1i*Egs*2*t);           % We multiply it by 2 the 
                                       % time-step. (1/N is for normal.)
                                       % Because for the greater Green 
                                       % function (eq.(7) Ref.[1]) we time
                                       % evolve only to tmax/2.
%}

%-------- tDMRG_star function

%{
function [M,EE,dw] = tDMRG_star(M,Hs,Nkeep,dt,tmax)
%{
< Description >


< Input >
 M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
       defines the chain length. The leg convention of M{n} is as follows:
    1      3   1      3         1        3
   ---M{1}---*---M{2}---* ... *---M{end}---
       |          |                 |
       ^2         ^2                ^2
 H : [cell] Hamiltonian. Each cell element H{n} describes the two-site
       interaction between site n and n+1. (Also has the information, given
       we are in star geometry, of the local energy). It should satisfy 
       numel(M) == numel(H) + 1.
       The leg convention of H{n} are as follows:
       2      4       [legs 1 and 2 are for site 
       |      |       legs 3 and 4 are for site n+1]
      [  H{n}  ]
       |      |
       1      3
 Nkeep : [integer] Maximum bond dimension.
 dt : [numeric] Real time step size. Each real-time evolution by step dt
       consists of texp(-dt/2*H) except exp(-dt*H) for the last step.
 tmax : [numeric] Maximum time range.
 Egs : [numeric] Ground state energy found using DMRG algorithm.

< Output >

%}   

tobj = tic2;

% Sanity Check
if length(M) ~= (length(Hs)+1)
    error("ERROR: it should be length(M) == (length(H)+1)");
end

N = numel(M);
Nstep = tmax/dt + 1;


EE = zeros(Nstep,(N-1));
dw = zeros(size(EE));

% SWAP gate
SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]; % matrix form in order to 
                                              % multiply it for exp(H).
                                        
% show information
fprintf("tDMRG: Real-time evolution - Star Gemotry\n");
fprintf(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', dt = ',sprintf('%i',dt),'\n']);

% generate the unitray operator exp(-it*H) for each two-site pairs

Fgate = cell(length(Hs)-1,1); % Forward gates (1st half-evolve then SWAP)
Rgate = cell(length(Hs)-1,1); % Reverse gates (1st SWAP then half-evolve)
Lgate = cell(1,1);            % Last gate (Full evolve only -  no SWAP's)


for it1 = (1:length(Hs))
    if ~isempty(Hs{it1})
        sdim = [size(M{it1},2),size(M{it1+1},2)];
        Htmp = permute(Hs{it1},[1 3 2 4]);
        Htmp = reshape(Htmp,[sdim(1)*sdim(2) sdim(1)*sdim(2)]);
        if it1 == length(Hs) % Last Hs acts with full time-step
            ttmp = dt; 
            eH = expm(-1i*ttmp*Htmp);
            Lgate{1} = reshape(eH, size(Hs{it1}));
        else
            ttmp = dt/2; % Half time step for all other steps
            eH = expm(-1i*ttmp*Htmp);
            Fgate{it1} = reshape(SWAP*eH, size(Hs{it1}));
            Rgate{it1} = reshape(eH*SWAP, size(Hs{it1})); 
        end
    end    
end

% bring into right-canonical form 
M = canonForm(M,0,'-k'); % one should not truncate zero singular values 
        % and the corresponding singular vectors, since it will decrease
        % the Hilbert space.

for it1 = (1:Nstep)
    % Left -> Right (Until end of Fgate)
    for it = (1:length(Fgate)) 
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,3,M{it+1},3,1);
        T = contract(Fgate{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        % SVD via svdTr
        [M{it},S,V,dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
        % normalize the singular values, to normalize the norm of MPS
        S = S/norm(S);
        % update M{it+1}
	    DS = diag(S);
        M{it+1} = contract(DS,2,2,V,3,1);
    end

    % Last step where we do NOT apply a SWAP gate
    T = contract(M{end-1},3,3,M{end},3,1);
    T = contract(Lgate{1},4,[3 4],T,4,[2 3],[3 1 2 4]);
    % SVD via svdTr
    [M{end-1},S,V,dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
    % normalize the singular values, to normalize the norm of MPS
    S = S/norm(S);
    % compute entanglement entropy of base 2. Be aware of zero
    % singular values that may produce Nan or Inf values!
    Evec = (S.^2).*log2(S.^2);
	for i = (1:length(Evec))
		if isnan(Evec(i)) || isinf(Evec(i))
			 Evec(i) = 0.0;
		end
	end
    EE(length(Hs)) = sum(-Evec);
    % update M{end}
	DS = diag(S);
    M{end} = contract(DS,2,2,V,3,1);
    M{end} = M{end}/norm(M{end}(:)); % to normalize the norm of MPS

    % Right -> Left (Using Rgate)
    for it = (length(Rgate):-1:1)
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,3,M{it+1},3,1);
        T = contract(Rgate{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        % SVD via svdTr
        [U,S,M{it+1},dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
        % normalize the singular values, to normalize the norm of MPS
        S = S/norm(S);
        % compute entanglement entropy of base 2. Be aware of zero
        % singular values that may produce Nan or Inf values!
        Evec = (S.^2).*log2(S.^2);
	    for i = (1:length(Evec))
		    if isnan(Evec(i)) || isinf(Evec(i))
			    Evec(i) = 0.0;
		    end
	    end
        EE(it1,it) = sum(-Evec); %2/log(2))*
        % update M{it}
        DS = diag(S);
        M{it} = contract(U,3,3,DS,2,1);
    end
    M{1} = M{1}/norm(M{1}(:)); % to normalize the norm of MPS (RIGHT)
    str = ['Sweep #',sprintf('%i/%i',[it1,Nstep])];
        disptime(str);
end
toc2(tobj,'-v');
end
%}
%{
'after right sweep'
XM = [];
for itN = (N:-1:1)
    Tf = permute(M{itN},[3 2 1]); % permute left<->right to use uLeft                                         
    Tr = permute(M{itN},[3 2 1]); % permute left<->right to use uLeft
    XM = updateLeft(XM,2,Tr,[],[],Tf);
end
XM
%}

%% Trials
%{
M = canonForm(M,Nsite,'-k');


dt = 0.1;
xt = 0:dt:20;
[U1,E1] = eig((Hchain+Hchain')/2);
Ud = U1';
Gexact= zeros(length(xt),1);
Ed = diag(E1);
Cn2 = zeros(Nsite,1);  % |C_n|^2
Crx = cell(Nsite,1);  % Creation MPOs.
Anx = cell(Nsite,1);  % Annihilation MPOs.
sm = squeeze(F);
sp = squeeze(F)';
sz = squeeze(Z);
for itk = 1:(Nsite)
    % First: build annihilation operator MPO.
    Anx{1} = zeros(1,2,2,2);
    Anx{1}(1,:,1,:) = Ud(itk,1)*sm;
    Anx{1}(1,:,2,:) = eye(2);
    for itj = 2:(Nsite-1)
        Anx{itj} = zeros(2,2,2,2);
        Anx{itj}(1,:,1,:) = sz;
        Anx{itj}(2,:,1,:) = Ud(itk,itj)*sm;
        Anx{itj}(2,:,2,:) = eye(2);
    end
    Anx{Nsite} = zeros(2,2,1,2);
    Anx{Nsite}(1,:,1,:) = sz;
    Anx{Nsite}(2,:,1,:) = Ud(itk,Nsite)*sm;
    % Second: build creation operator MPO.
    Crx{1} = zeros(1,2,2,2);
    Crx{1}(1,:,1,:) = U1(1,itk)*sp;
    Crx{1}(1,:,2,:) = eye(2);
    for iti = 2:(Nsite-1)
        Crx{iti} = zeros(2,2,2,2);
        Crx{iti}(1,:,1,:) = sz;
        Crx{iti}(2,:,1,:) = U1(iti,itk)*sp;
        Crx{iti}(2,:,2,:) = eye(2);
    end
    Crx{Nsite} = zeros(2,2,1,2);
    Crx{Nsite}(1,:,1,:) = sz;
    Crx{Nsite}(2,:,1,:) = U1(Nsite,itk)*sp;
    % Apply the number operator to Psi0.
        % % contract the left half.
    cl = contract(Anx{1},4,4,M{1},3,2);
    cl = contract(Crx{1},4,4,cl,5,2);
    cl = contract(conj(M{1}),3,2,cl,7,2);
    cl = permute(cl,[2,4,6,8,1,3,5,7]);
    for j = 2:(Nsite)
        cl = contract(cl,4,4,M{j},3,1);
        cl = contract(Anx{j},4,[1,4],cl,5,[3,4]);
        cl = contract(Crx{j},4,[1,4],cl,5,[4,1]);
        cl = contract(conj(M{j}),3,[1,2],cl,5,[4,1]);
    end
    Cn2(itk) = cl;
    % % combine left and right halves.
    %Cn2(itk) = contract(cl,4,[1,2,3,4],cr,4,[1,2,3,4]);
end
anorm = sum(Cn2);
for itT = 1:numel(xt)
    Gexact(itT) = sum(Cn2.*exp(-1i*(Ed-E0)*xt(itT)))/anorm;
end
figure
plot(xt,real(Gexact))
grid on
%}

% Construction of the Hs (2site) for the star geometry

% -#######################################################################-

%% FUNCTIONS
%{ 
From here the code contains the functions needed to to the plot, except the
ones that were already in the code repository. Which I assume that are
already a given.
%}

% -#######################################################################-
% -#######################################################################-

%% doLinDisc

function [ff,gg,estarimp,vstarimp] = doLinDisc(ozin,RhoV2in,N)
%{
< Description >


< Input >


< Output >

%}

deltaE = ozin(2)-ozin(1);
estar = zeros(1,N);
vstar = zeros(1,N);

% N equidistant intervals Ik of size ∆e. 
% e_k = min{Ik} + ∆e/2
for i=(1:N)
    estar(i) = min(ozin(i),ozin(i+1)) + deltaE/2;
end

for i=(1:N)
    a = [ozin(i) ozin(i+1)];
    b = [RhoV2in(i) RhoV2in(i+1)];
    vstar(i) = sqrt(trapz(a,b));
end

%{
estarimp = [0 estar(N:-1:(N+3)/2) estar(1:(N-1)/2) estar((N+1)/2)];
vstarimp = [0 vstar(N:-1:(N+3)/2) vstar(1:(N-1)/2) vstar((N+1)/2)];
%}

%
estarimp = [0 estar(1:N)];
vstarimp = [0 vstar(1:N)];
%}
Hstar = diag(estarimp);
Hstar(:,1) = vstarimp;
Hstar(1,:) = vstarimp;


U = zeros(size(Hstar,1),1);
U(1,1) = 1;

% Lanczos tridiagonalization
ff = zeros(N,1); % hopping amplutudes; corresponds to the super diagonal 
                  % (diag(..,+1) or diag(..,-1)) in the tridiagonal matrix
gg = zeros(N+1,1); % hopping amplutudes; corresponds to the main diagonal 
                  % (diag(..)) in the tridiagonal matrix


for itN = (1:N)
    v = Hstar*U(:,itN);
    v = v-U*(U'*v);
    v = v-U*(U'*v); % twice for numerical reason
    ff(itN) = norm(v);

    if (itN < N) && (ff(itN) > 0)
        U(:,itN+1) = v/ff(itN);
        gg(itN) = U(:,itN+1)'*Hstar*U(:,itN+1);
    elseif (itN < N) && (ff(itN) <= 0)
        ff = ff(1:itN);
        gg = gg(1:itN);
        break
    end
end

figure
set(gcf,'units','points');
set(gcf,'position',[0 0 450 600]);
set(gcf,'paperunits', 'points');
set(gcf,'papersize', [450 600]);
set(gcf,'paperpositionmode', 'manual');
set(gcf,'paperposition', [0 0 450 600]);
set(gca,'units','points');
set(gca,'innerposition',[50,50,350,500]);
plot((1:N),ff,'+','linewidth',2,'MarkerSize',10);
hold on
plot((1:N),6*vstar,'x','linewidth',2,'MarkerSize',10);
hold off
xlabel('$\ell$','Interpreter','latex','fontsize',18);
xlim([0 N+0.5])
set(gca,'FontSize',18);
set(gca,'yminortick','on');
set(gca,'ticklength',[0.04 0.1]);
legend('$t_{chain}$',...
    '$6\cdot t_{star}$', ...
    'Interpreter','latex','fontsize',18,'location','south');

figure
set(gcf,'units','points');
set(gcf,'position',[0 0 450 600]);
set(gcf,'paperunits', 'points');
set(gcf,'papersize', [450 600]);
set(gcf,'paperpositionmode', 'manual');
set(gcf,'paperposition', [0 0 450 600]);
set(gca,'units','points');
set(gca,'innerposition',[50,50,350,500]);
plot((1:N+1),gg,'+','linewidth',2,'MarkerSize',10);
hold on
plot((1:N),estar,'x','linewidth',2,'MarkerSize',10);
hold off
xlabel('$\ell$','Interpreter','latex','fontsize',18);
xlim([0 N+0.5])
set(gca,'FontSize',18);
set(gca,'yminortick','on');
set(gca,'ticklength',[0.04 0.1]);
legend('$\varepsilon_{chain}$',...
    '$\varepsilon_{star}$', ...
    'Interpreter','latex','fontsize',18,'location','south');
end

% -#######################################################################-
% -#######################################################################-

%% tDMRG

function [Mf,Mr,EE,dw,Green] = tDMRG(M,Hs,Nkeep,dt,tmax)
%{
 < Description >
 [ts,M,EE,dw] = tDMRG (M,Hs,Nkeep,dt,tmax)
 Real time evolution by using time-dependent DMRG (tDMRG) for simulating
 real-time evolution of matrix product state (MPS). The 1D chain system is 
 described by the Hamiltonian H. Here we use the 2nd order Trotter 
 step exp(-dt/2*Hodd) * exp(-dt*Heven) * exp(-dt/2*Hodd). After acting each 
 exp(-t*H) type operator, the bonds are truncated such that the largest 
 singular values are kept. Also this function directly computes the overlap
 of the forward and reverse time evolutions, needed to then calculate the 
 greater Green's function.
 < Input >
 M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
       defines the chain length. The leg convention of M{n} is as follows:
    1      3   1      3         1        3
   ---M{1}---*---M{2}---* ... *---M{end}---
       |          |                 |
       ^2         ^2                ^2
 Hs : [cell] Hamiltonian. Each cell element H{n} describes the two-site
       interaction between site n and n+1. Thus, H(1:2:end) acts on odd
       bonds; H(2:2:end) on even bonds. It should satisfy numel(M) ==
       numel(H) + 1.
       The leg convention of H{n} are as follows:
       2      4       [legs 1 and 2 are for site n;
       |      |       legs 3 and 4 are for site n+1]
      [  H{n}  ]
       |      |
       1      3
 Nkeep : [integer] Maximum bond dimension.
 dt : [numeric] Real time step size. Each real-time evolution by step dt
       consists of three Trotter steps, exp(-dt/2*Hodd) * exp(-dt*Heven) *
       exp(-dt/2*Hodd).
 tmax : [numeric] Maximum time range.
 Egs : [numeric] Ground state energy found using DMRG algorithm.
 < Output >
 Mf(Mr) : [cell] The final MPS's after real-time evolution for the two 
          evolutions needen in the greater Green function calculation.
 EE : [matrix] EE(m,n) indicates the entanglement entropy (with base 2) of
       the MPS with respect to the bipartition at the bond between the
       sites n and n+1, after acting the m-th Trotter step. (log2 -> ebits)
 dw : [matrix] Discarded weights (i.e., the sum of the squares of the
       discarded singular values) after each Trotter step. dw(m,n)
       corresponds to the same bond and Trotter step associated with
       EE(m,n).
 Green: [vector] value of the greater Green function for each time step.

%}   

tobj = tic2;

%
% Sanity Check
if length(M) ~= (length(Hs)+1)
    error("ERROR: it should be length(M) == (length(H)+1)\n");
end
%}
N = numel(M);
Nstep = tmax/dt;

EE = zeros(Nstep,N-1);
dw = zeros(size(EE));
Green = zeros(Nstep,1);

% show information
fprintf("\n tDMRG: Real-time evolution - Chain Gemotry\n");
fprintf(['# of sites = ',sprintf('%i',numel(M)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', dt = ',sprintf('%i',dt),'\n']);

% generate the unitray operator exp(-it*H) for each two-site pairs
expH = cell(N-1,1);

% generate the unitray operator exp(+it*H) for each two-site pairs
expHr = cell(N-1,1);

for it1 = (1:length(Hs))
    if ~isempty(Hs{it1})
        sdim = [size(M{it1},2),size(M{it1+1},2)];
        Htmp = permute(Hs{it1},[1 3 2 4]);
        Htmp = reshape(Htmp,[sdim(1)*sdim(2) sdim(1)*sdim(2)]);
        if mod(it1,2) % True when odd
            ttmp = dt/2; % half time step for odd bonds, as the time 
            % evolution for odd bonds will happen twice (TS 2nd order)
        else
            ttmp = dt;
        end
        % Forward Time Evolution gates
        eH = expm(-1i*ttmp*Htmp);
        expH{it1} = reshape(eH, size(Hs{it1}));
        % Reverse Time Evolution gates
	    eHr = expm(1i*ttmp*Htmp);
        expHr{it1} = reshape(eHr, size(Hs{it1}));
    end    
end

% bring into right-canonical form 
M = canonForm(M,0,'-k'); % one should not truncate zero singular values 
        % and the corresponding singular vectors, since it will decrease
        % the Hilbert space.

Mf = M;
Mr = M;
for it1 = (1:3*Nstep)
% Here we use the 2nd order Trotter step exp(-dt/2*Hodd) * exp(-dt*Heven) *
% exp(-dt/2*Hodd). That is, for the case mod(it1,3) == 2, we act the
% unitary on even bonds. Otherwise, on odd bonds.
% variables ending with r refer to the reverse time evolution needed for
% the calculation of the greater Green function.
    expHtmp = cell(N-1,1);
    expHtmpr = cell(N-1,1);
    if mod(it1,3) == 2 	% even bonds
        for i = (2:2:numel(expH))
            expHtmp{i} = expH{i};
            expHtmpr{i} = expHr{i};
        end
    else 			    % odd bonds
        for i = (1:2:numel(expH))
            expHtmp{i} = expH{i};
            expHtmpr{i} = expHr{i};
        end
    end   
    % call local function tDMRG_1sweep
    [Mf,EE1,dw1] = tDMRG_1sweep(Mf,expHtmp,Nkeep,mod(it1,2));
    [Mr,~,~] = tDMRG_1sweep(Mr,expHtmpr,Nkeep,mod(it1,2));   
    
    % update the rows of entanglement entropy (EE) and discarded weights 
    % We do it every 3 iteration (full ∆t evolution).
    if mod(it1,3) == 0
        EE(it1/3,:) = EE1;
        dw(it1/3,:) = dw1;
        Green(it1/3) = tDMRG_expVal(Mf,Mr,mod(it1,2));
        
        %{
        MM = []; % contraction of bra/ket tensors
        if mod(it1,2) ~= 1 	% right-normalized
            for itN = (N:-1:1)
                Tf = permute(Mf{itN},[3 2 1]); % permute left<->right                               
                Tr = permute(Mr{itN},[3 2 1]); % permute left<->right 
                MM = updateLeft(MM,2,Tf,[],[],Tr);
            end
        else 	% left-normalized
            for itN = (1:N)
                 MM = updateLeft(MM,2,Mf{itN},[],[],Mr{itN});
            end
            Green(it1/3) = MM*exp(1i*2*dt*Egs);
        end
        %}
    end
    % display informaiton of the sweep
    if mod(it1,3) == 0
        str = ['Sweep #',sprintf('%i/%i',[it1/3,Nstep])];
        disptime(str);
    end
end

       
toc2(tobj,'-v');
end

% -------------------------------------------------------------------------

function [M,EE,dw] = tDMRG_1sweep(M,expH,Nkeep,isright)
%{
 Apply exp(-it*H), which is an array of two-site gates acting on either
 even or odd bonds, and then truncate bonds by using SVD. After applying
 this function, left-canonical state becomes right-canonical, and vice
 versa.
 This is what we called in class bond-by-bond compression.
 < Input >
 M : [cell] Input MPS.
 expH : [cell] exp(-i*H*T) unitary operators for each bond. The length
       should satisfy numel(expH) == numel(M)-1. And the every first (or
       second) elements should be empty, since we act either even or odd
       bonds at once.
 Nkeep : [numeric] Maximum bond dimension.
 isright : [logical] If true, we take left-to-right sweep. Otherwise, take
       right-to-left sweep.
 
 < Output >
 M : [cell] MPS after applying exp(-it*H) and truncating bonds.
 EE : [numeric vector] Entanglement entropy at each bond.
 dw : [numeric vector] Discarded weights when truncating the bond
       dimensions.
%}

N = length(M);
EE = zeros(1,N-1);
dw = zeros(1,N-1);

if isright == 1 % left -> right
    for it = (1:N-1)
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,3,M{it+1},3,1);
        if ~isempty(expH{it})
            T = contract(expH{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        end
        % SVD via svdTr
        [M{it},S,V,dw(it)] = svdTr(T,4,[1,2],Nkeep,-1);
        % normalize the singular values, to normalize the norm of MPS
        S = S/norm(S);
        % compute entanglement entropy of base 2. We take care of 0
        % singular values that may produce Nan or Inf values!
        Evec = (S.^2).*log2(S.^2);
	for i = (1:length(Evec))
		if isnan(Evec(i)) || isinf(Evec(i))
			Evec(i) = 0.0;
		end
	end
        EE(it) = sum(-Evec);
        % update M{it+1}
	    DS = diag(S);
        M{it+1} = contract(DS,2,2,V,3,1);
    end
    M{end} = M{end}/norm(M{end}(:)); % to normalize the norm of MPS

else % right -> left
    for it = (N-1:-1:1)
        % contract M{it} and M{it+1} with expH{it}
        T = contract(M{it},3,3,M{it+1},3,1);
        if ~isempty(expH{it})

            T = contract(expH{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
        end
        % SVD via svdTr
        [U,S,M{it+1},dw(it)] = svdTr(T,4,[1,2],Nkeep,-1);
        % normalize the singular values, to normalize the norm of MPS
        S = S/norm(S);
        % compute entanglement entropy of base 2. We take care of 0
        % singular values that may produce Nan or Inf values!
        Evec = (S.^2).*log2(S.^2);
	for i = (1:length(Evec))
		if isnan(Evec(i)) || isinf(Evec(i))
			Evec(i) = 0.0;
		end
	end
        EE(it) = sum(-Evec); %2/log(2))*
        % update M{it}
        DS = diag(S);
        M{it} = contract(U,3,3,DS,2,1);
    end
    M{1} = M{1}/norm(M{1}(:)); % to normalize the norm of MPS
end

end

% -------------------------------------------------------------------------

function Green = tDMRG_expVal(Mf,Mr,isleft)
%{
 Computes the overlap between the forward and reverse time evolutions.
 < Input >
 Mf(Mr) : [cell] Input MPS's.
 Egs : [numerical] Ground State Energy found by 2s DMRG.
 isleft : [logical] If true (==1), it means that the MPS M is in 
       left-canonical form. Otherwise, right-canonical form.
 dt : [numerical] Time-Step
 < Output >
 Ovals : [vector] Green function at the desired time step.
%}

N = length(Mf);

% Sanity Check
if N ~= length(Mr)
    error('Something went wrong!')
end

%Ovals = zeros(N,1);

MM = []; % contraction of bra/ket tensors
if isleft ~= 1 	% right-normalized
    for itN = (N:-1:1)
        Tf = permute(Mf{itN},[3 2 1]); % permute left<->right to use uLeft                                         
        Tr = permute(Mr{itN},[3 2 1]); % permute left<->right to use uLeft
        MM = updateLeft(MM,2,Tr,[],[],Tf);
    end
else 	% left-normalized
    for itN = (1:N)
        MM = updateLeft(MM,2,Mr{itN},[],[],Mf{itN});
    end
end

Green = MM;

end

% -#######################################################################-
% -#######################################################################-

%% tDMRG_star

function [Mf,Mr,EE,dw,Green] = tDMRG_star(M,Hs,Nkeep,dt,tmax)
%{
< Description >


< Input >
 M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
       defines the chain length. The leg convention of M{n} is as follows:
    1      3   1      3         1        3
   ---M{1}---*---M{2}---* ... *---M{end}---
       |          |                 |
       ^2         ^2                ^2
 H : [cell] Hamiltonian. Each cell element H{n} describes the two-site
       interaction between site n and n+1. (Also has the information, given
       we are in star geometry, of the local energy). It should satisfy 
       numel(M) == numel(H) + 1.
       The leg convention of H{n} are as follows:
       2      4       [legs 1 and 2 are for site 
       |      |       legs 3 and 4 are for site n+1]
      [  H{n}  ]
       |      |
       1      3
 Nkeep : [integer] Maximum bond dimension.
 dt : [numeric] Real time step size. Each real-time evolution by step dt
       consists of texp(-dt/2*H) except exp(-dt*H) for the last step.
 tmax : [numeric] Maximum time range.
 Egs : [numeric] Ground state energy found using DMRG algorithm.

< Output >

%}   

tobj = tic2;

% Sanity Check
if length(M) ~= (length(Hs)+1)
    error("ERROR: it should be length(M) == (length(H)+1)");
end

N = numel(M);
Nstep = tmax/dt;


EE = zeros(Nstep,N-1);
dw = zeros(size(EE));
Green = zeros(Nstep,1);

% SWAP gate
SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]; % matrix form in order to 
                                              % multiply it for exp(H).
                                        
% show information
fprintf("\n tDMRG: Real-time evolution - Star Gemotry\n");
fprintf(['# of sites = ',sprintf('%i',numel(M)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', dt = ',sprintf('%i',dt),'\n']);

% generate the unitray operator exp(-it*H) for each two-site pairs

Fgate = cell(length(Hs)-1,1); % Forward gates (1st half-evolve then SWAP)
Rgate = cell(length(Hs)-1,1); % Reverse gates (1st SWAP then half-evolve)
Lgate = cell(1,1);            % Last gate (Full evolve only -  no SWAP's)

% Reverse evolution gates
Fgater = cell(length(Hs)-1,1); % Forward gates (1st half-evolve then SWAP)
Rgater = cell(length(Hs)-1,1); % Reverse gates (1st SWAP then half-evolve)
Lgater = cell(1,1);            % Last gate (Full evolve only -  no SWAP's)


for it1 = (1:length(Hs))
    if ~isempty(Hs{it1})
        sdim = [size(M{it1},2),size(M{it1+1},2)];
        Htmp = permute(Hs{it1},[1 3 2 4]);
        Htmp = reshape(Htmp,[sdim(1)*sdim(2) sdim(1)*sdim(2)]);
        if it1 == length(Hs) % Last Hs acts with full time-step
            ttmp = dt; 
            eH = expm(-1i*ttmp*Htmp);
            Lgate{1} = reshape(eH, size(Hs{it1}));
            eHr = expm(+1i*ttmp*Htmp);
            Lgater{1} = reshape(eHr,size(Hs{it1}));
        else
            ttmp = dt/2; % Half time step for all other steps
            eH = expm(-1i*ttmp*Htmp);
            Fgate{it1} = reshape(SWAP*eH, size(Hs{it1}));
            Rgate{it1} = reshape(eH*SWAP, size(Hs{it1})); 
            eHr = expm(+1i*ttmp*Htmp);
            Fgater{it1} = reshape(SWAP*eHr, size(Hs{it1}));
            Rgater{it1} = reshape(eHr*SWAP, size(Hs{it1})); 
        end
    end    
end

% bring into right-canonical form 
M = canonForm(M,0,'-k'); % one should not truncate zero singular values 
        % and the corresponding singular vectors, since it will decrease
        % the Hilbert space.


Mf = M;
Mr = M;

for it1 = (1:Nstep)
    % call local function S_sweep
    %[M,EE,dw] = S_sweep(M,Fgate,Rgate,Lgate,Nkeep)
    [Mf,EE1,dw1] = S_sweep(Mf,Fgate,Rgate,Lgate,Nkeep);
    [Mr,~,~] = S_sweep(Mr,Fgater,Rgater,Lgater,Nkeep);

    EE(it1,:) = EE1;
    dw(it1,:) = dw1;
    Green(it1) = tDMRG_expValS(Mf,Mr);

    str = ['Sweep #',sprintf('%i/%i',[it1,Nstep])];
    disptime(str);
end

toc2(tobj,'-v');
end

% -------------------------------------------------------------------------

function [M,EE,dw] = S_sweep(M,Fgate,Rgate,Lgate,Nkeep)
N = numel(M);
 % Left -> Right (Until end of Fgate)
for it = (1:length(Fgate)) 
    % contract M{it} and M{it+1} with expH{it}
    T = contract(M{it},3,3,M{it+1},3,1);
    T = contract(Fgate{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
    % SVD via svdTr
    [M{it},S,V,dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
    % normalize the singular values, to normalize the norm of MPS
    S = S/norm(S);
    % update M{it+1}
	DS = diag(S);
    M{it+1} = contract(DS,2,2,V,3,1);
end

% Last step where we do NOT apply a SWAP gate
T = contract(M{N-1},3,3,M{N},3,1);
T = contract(Lgate{1},4,[3 4],T,4,[2 3],[3 1 2 4]);
% SVD via svdTr
[M{end-1},S,V,dw(N-1)] = svdTr(T,4,[1,2],Nkeep,[]);
% normalize the singular values, to normalize the norm of MPS
S = S/norm(S);
% compute entanglement entropy of base 2. Be aware of zero
% singular values that may produce Nan or Inf values!
Evec = (S.^2).*log2(S.^2);
for i = (1:length(Evec))
    if isnan(Evec(i)) || isinf(Evec(i))
	    Evec(i) = 0.0;
    end
end
EE(N-1) = sum(-Evec);
% update M{end}
DS = diag(S);
M{end} = contract(DS,2,2,V,3,1);
M{end} = M{end}/norm(M{end}(:)); % to normalize the norm of MPS

% Right -> Left (Using Rgate)
for it = (length(Rgate):-1:1)
    % contract M{it} and M{it+1} with expH{it}
    T = contract(M{it},3,3,M{it+1},3,1);
    T = contract(Rgate{it},4,[3 4],T,4,[2 3],[3 1 2 4]);
    % SVD via svdTr
    [U,S,M{it+1},dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
    % normalize the singular values, to normalize the norm of MPS
    S = S/norm(S);
    % compute entanglement entropy of base 2. Be aware of zero
    % singular values that may produce Nan or Inf values!
    Evec = (S.^2).*log2(S.^2);
	for i = (1:length(Evec))
		if isnan(Evec(i)) || isinf(Evec(i))
			Evec(i) = 0.0;
		end
	end
    EE(it) = sum(-Evec); %2/log(2))*
    % update M{it}
    DS = diag(S);
    M{it} = contract(U,3,3,DS,2,1);
end
M{1} = M{1}/norm(M{1}(:)); % to normalize the norm of MPS (RIGHT)

end

% -------------------------------------------------------------------------

function Green = tDMRG_expValS(Mf,Mr)
%{
 Computes the overlap between the forward and reverse time evolutions.
 < Input >
 Mf(Mr) : [cell] Input MPS's.
 < Output >
 Ovals : [vector] Overlap of MPS's at a given time step.
%}

N = length(Mf);

% Sanity Check
if N ~= length(Mr)
    error('Something went wrong!')
end

Mf = canonForm(Mf,0,'-k');
Mr = canonForm(Mr,0,'-k');

MM = []; % contraction of bra/ket tensors
for itN = (N:-1:1)
    Tf = permute(Mf{itN},[3 2 1]); % permute left<->right to use uLeft                                         
    Tr = permute(Mr{itN},[3 2 1]); % permute left<->right to use uLeft
    MM = updateLeft(MM,2,Tr,[],[],Tf);
end

Green = MM;

end

