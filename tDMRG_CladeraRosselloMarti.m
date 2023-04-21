% FRAME 2: 2. tDMRG for single-impurity Anderson model (SIAM)
% MSc Quantum Science and Technology - Tensor Networks
% Martí Cladera Rosselló



%% IMPORTANT
% The code is created so that a window pups-up to let you choose if we want
% to study the system in the case of no interaction (U=0) or interaction
% (U=1). Also another menu lets you choose if you want to see the plots 
% of the linear discretization of the bath (doLinDisc.m). Finally another 
% pop-up menu lets you choose if you want to see the plots from the time 
% evolution and entropy or not.

clear
clc

N = 9; % Number of Bath Sites

% I created the code so that when having U = 0 we exploit the fact that the
% system becomes spin independent (There's no interaction between spin up
% and spin down species in our MPS description), and hence, the
% computations become cheaper.

choice = menu(['Do you want to study the Non-Interacting [U=0]' ...
    ' or Interacting [U=1] case?'],'U = 1','U = 0');
if choice==2 || choice==0
    U = 0;
else
    U = 1;
end

epsd = -U/2; % On-site Impurtiy energy
Nkeep = 250; % Maximum Bond Dimension
dt = 0.1; % 2*∆t
tmax = 25; 
Tstep = 0:dt:tmax;

% Initialize Local space of "spinless fermions". The fact that we consider
% them spinless is what doubles the size of our MPS,
[F,Z,I] = getLocalSpace('Fermion');

if U == 0
    Nsite = (N+1);
else
    Nsite = 2*(N+1);
end
fprintf("\n 2. tDMRG for single-impurity Anderson model (SIAM)\n" + ...
    "Martí Cladera Rosselló \n");
fprintf(['\n Parameters: U = ',sprintf('%i',U), ...
    ', epsd = ',sprintf('%.2g',epsd),...
    ', Effective Nsite = ',sprintf('%i',Nsite),'\n\n']);


%% Discretization of the bath and Hamiltonian creation

ozin = linspace(-1,1,N+1);
[ff,gg,estarimp,vstarimp] = doLinDisc(ozin,N);
% ff(gg) are the hopping(on-site) parameters in chain geometry for U=0

estar = estarimp(2:end); % On-site bath energies - Star Geometry (U=0)
vstar = vstarimp(2:end); % Hopping parameters - Star Geometry (U=0)
ff_sp = [flipud(ff);0;ff]; % Hopping parameters - Chain Geometry (U=1)
gg_sp = [flipud(gg);gg]; % On-site bath energies - Chain Geometry (U=1)
estarimp_sp = [fliplr(estar) epsd epsd estar]; % On-site bath energies
                                           %  Star Geometry (U=1)
vstarimp_sp = [fliplr(vstarimp) vstarimp]; % Hopping parameters  
                                           %  Star Geometry (U=1)


% We construct the Chain and Star Hamiltonians for the non-interacting case 
% to then compute the exact Green Function.

if U == 0
    % Construction of H_star
    Hstar = diag(estarimp);
    Hstar(:,1) = vstarimp;
    Hstar(1,:) = vstarimp;

    % Construction of H_chain
    Hchain = diag(ff,1);
    Hchain = Hchain +Hchain' + diag(gg);

end

%% Ground State search - DMRG

% Let us construct the MPO, notice how the Hchain has the form of 
% SUM_{nn}(t_{ij}*C^*_{i}*C_{j} + h.c) HOWEVER we don't have to take 
% into account the spin degree of freedom due to the particle-hole symmetry
% so the Hilbert space is two dimensional (empty or occupied).
% Notice how we number of bath sites : 2N -> 2(N+1) counting the impurity.
% Remark 1: We'll be following the MPS description in Ref.[1].
% Remark 2: When U = 0 we can neglect the differenciation between 
% spin-up and spin-down sites, since the Hamiltonian becomes spin
% independent. We'll be doing the calculations for only one spin when U=0.
%{
u: spin-up
d: spin-down
u-I: spin-up Impurity
d-I: spin-down Impurity
Ref[1]
 
         BATH-u                                           BATH-d
    --------^-------------                       -----------^-----------
Site:  N+1        N            1           1            N        N+1

    ---(u)---*---(u)---...---(u-I)---*---(d-I)---...---(d)---*---(d)---
        |         |            |           |            |         |
        |         |            |           |            |         |
Comp    
Site:   1         2           N+1         N+2         2(N+1)    2(N+2)

%}

if U ~= 0
    % As stated in the paper by Daniel Bauernfeind et.al., we'll be using a
    % local Hilbert of dimension 2 (empty |0> or occupied |1>)
    % Here we take into account the system for U~=0, where we cannot treat 
    % the system as an N+1 chain, but a 2(N+1) chain has to be considered.

     W = cell(Nsite,1); % MPO W
     % First site
     W{1} = zeros(1,2,4,2); % ordering: left bottom right top
     W{1}(1,:,2,:) = ff_sp(1)*squeeze(F)'; % t_l*creation
     W{1}(1,:,3,:) = ff_sp(1)'*squeeze(F); % t_l*annihilation
     W{1}(1,:,4,:) = I;
        % Last site
     W{Nsite} = zeros(4,2,1,2); % ordering: left bottom right top
     W{Nsite}(1,:,1,:) = I;
     W{Nsite}(2,:,1,:) = squeeze(F); % annihilation
     W{Nsite}(3,:,1,:) = squeeze(F)'; % creation
     % Other sites
     for itL2 = (2:N)
            W{itL2} = zeros(4,2,4,2); % ordering: left bottom right top
            W{itL2}(1,:,1,:) = I;
            W{itL2}(2,:,1,:) = squeeze(F); % annihilation
            W{itL2}(3,:,1,:) = squeeze(F)'; % creation
            W{itL2}(4,:,2,:) = ff_sp(itL2)*squeeze(F)'; % t_l*creation
            W{itL2}(4,:,3,:) = ff_sp(itL2)'*squeeze(F); % t_l*annihilation
            W{itL2}(4,:,4,:) = I;
     end
     W{N+1} = zeros(4,2,4,2);
     W{N+1}(1,:,1,:) = I;
     W{N+1}(2,:,1,:) = squeeze(F); % annihilation
     W{N+1}(3,:,1,:) = squeeze(F)'; % creation
     W{N+1}(4,:,1,:) = epsd*squeeze(F)'*squeeze(F); % epsd*number
     W{N+1}(4,:,2,:) = U*squeeze(F)'*squeeze(F); % U*number
     W{N+1}(4,:,4,:) = I;
    
     W{N+2} = zeros(4,2,4,2);
     W{N+2}(1,:,1,:) = I;
     W{N+2}(2,:,1,:) = squeeze(F)'*squeeze(F); % number
     W{N+2}(4,:,1,:) = epsd*squeeze(F)'*squeeze(F); % epsd*number
     W{N+2}(4,:,2,:) = ff_sp(N+2)*squeeze(F)'; % t_I*creation
     W{N+2}(4,:,3,:) = ff_sp(N+2)*squeeze(F); % t_I*annihilation
     W{N+2}(4,:,4,:) = I;

     for itL2 = ((N+3):(Nsite-1))
            W{itL2} = zeros(4,2,4,2); % ordering: left bottom right top
            W{itL2}(1,:,1,:) = I;
            W{itL2}(2,:,1,:) = squeeze(F); % annihilation
            W{itL2}(3,:,1,:) = squeeze(F)'; % creation
            W{itL2}(4,:,2,:) = ff_sp(itL2)*squeeze(F)'; % t_l*creation
            W{itL2}(4,:,3,:) = ff_sp(itL2)'*squeeze(F); % t_l*annihilation
            W{itL2}(4,:,4,:) = I;
     end

else
    % MPO for Hff: START %
    % As stated in the paper by Daniel Bauernfeind et.al., we'll be using a
    % local Hilbert of dimension 2 (empty |0> or occupied |1>)
    % If  U = 0 the MPS description followed in Ref.[1] becomes identical 
    % in both spin-up and spin-down sites, consequently, we can do it for 
    % just one of the spins descriptions, the Greater Green's function for
    % the other spin will be exactly the same, due to symmetry. Also 
    % notice how this reduces considerably the computational resources due 
    % to the the fact that we are just evolving half of the original system

     W = cell(Nsite,1); % MPO W
     % First site
     W{1} = zeros(1,2,4,2); % ordering: left bottom right top
     W{1}(1,:,2,:) = ff(1)*squeeze(F)'; % t_l*creation
     W{1}(1,:,3,:) = ff(1)*squeeze(F); % t_l*annihilation
     W{1}(1,:,4,:) = I;
        % Last site
     W{Nsite} = zeros(4,2,1,2); % ordering: left bottom right top
     W{Nsite}(1,:,1,:) = I;
     W{Nsite}(2,:,1,:) = squeeze(F); % annihilation
     W{Nsite}(3,:,1,:) = squeeze(F)'; % creation
     % Other sites
     for itL2 = (2:(Nsite-1))
            W{itL2} = zeros(4,2,4,2); % ordering: left bottom right top
            W{itL2}(1,:,1,:) = I;
            W{itL2}(2,:,1,:) = squeeze(F); % annihilation
            W{itL2}(3,:,1,:) = squeeze(F)'; % creation
            W{itL2}(4,:,2,:) = ff(itL2)*squeeze(F)'; % t_l*creation
            W{itL2}(4,:,3,:) = ff(itL2)*squeeze(F); % t_l*annihilation
            W{itL2}(4,:,4,:) = I;
     end
end

MPS = cell(Nsite,1);
MPS{1} = rand(1,2,Nkeep);
MPS(2:Nsite-1) = {rand(Nkeep,2,Nkeep)};
MPS{Nsite} = rand(Nkeep,2,1);


if U == 0
    [M0,E0,~] = DMRG_2site(W,MPS,sqrt(2),Nkeep,100,'Econv',1e-12);
    fprintf(['\n  Non-Interacting Ground State Energy (U = 0) = ' ...
        ,sprintf('%.7g',E0),'\n'])
else
    [M0ref,E0ref,~] = DMRG_2site(W,MPS,sqrt(2),Nkeep,100,'Econv',1e-14);
    [M0,E0,~] = DMRG_2site(W,MPS,sqrt(2),Nkeep/2,100,'Econv',1e-12);
    fprintf(['\n Interacting Ground State Energy (U = 1) = ', ...
        sprintf('%.7g',E0ref),'\n'])
end

%% Creating the initial state

% The initial state is applying the creation operator on the impurity site
% of the ground state.

Fdag = squeeze(F)';
M = cell(Nsite,1);

if U == 0
    % Let us construct the intitial MPO applying C^dad_I only to the first 
    % site corresponding to the impurity, and the corresponding Z to take 
    % into account the fermionic signs:

    M{1} = contract(Fdag,2,2,M0{1},3,2,[2 1 3]);
    for i=(2:(numel(M0)))
        M{i} = contract(Z,2,2,M0{i},3,2,[2 1 3]); %M0{i};
    end
else
    % Due to the MPS description used throughout the project and in the
    % reference, I was unsure on how to apply the creation operators. I was
    % doubting between the usual ordering of the MPS {1...2(N+1)}, which is
    % the one I ended up using, or if I had to apply it as
    % {(N+1)...1,1...(N+1)}. In the latter, I would have to had applied Z
    % gates to all the non-impurity sites (done in the comment below).

    Mref = M0;
    Mref{N+1} = contract(Fdag,2,2,M0ref{N+1},3,2,[2 1 3]);
    Mref{N+2} = contract(Z*Fdag,2,2,M0ref{N+2},3,2,[2 1 3]);

    M = M0;
    M{N+1} = contract(Fdag,2,2,M0{N+1},3,2,[2 1 3]);
    M{N+2} = contract(Z*Fdag,2,2,M0{N+2},3,2,[2 1 3]);

    %{
    % Reference
    Mref{N+1} = contract(Fdag,2,2,M0ref{N+1},3,2,[2 1 3]);
    Mref{N+2} = contract(Fdag,2,2,M0ref{N+2},3,2,[2 1 3]);
    for i=(1:N)
        Mref{i} = contract(Z,2,2,M0ref{i},3,2,[2 1 3]);
    end
    for i=(N+3:length(M0ref))
        Mref{i} = contract(Z,2,2,M0ref{i},3,2,[2 1 3]);
    end

    % Less precision
    M{N+1} = contract(Fdag,2,2,M0{N+1},3,2,[2 1 3]);
    M{N+2} = contract(Fdag,2,2,M0{N+2},3,2,[2 1 3]);
    for i=(1:N)
        M{i} = contract(Z,2,2,M0{i},3,2,[2 1 3]);
    end
    for i=(N+3:length(M0))
        M{i} = contract(Z,2,2,M0{i},3,2,[2 1 3]);
    end
    %}
end


%% Creating the two-site Hamiltonians

% Construct the Hamiltonian as two site operators for Chain geometry
Hc = cell(Nsite-1,1);

n_I = [0 0; 0 1]; % Number operator 1 if occupied (|1><1|) 0 otherwise.

if U == 0
    for i =(1:length(ff))
        % gg(i) = 0, no on-site bath term in any site of the chain
        Hc{i} = ff(i)*contract(F,3,2,permute(conj(F),[3 2 1]),3,2);
        Hc{i} = Hc{i} + ff(i)*contract(permute(conj(F), ...
            [3 2 1]),3,2,F,3,2);
    end
else 
    Tensor_N1 =reshape(epsd*(kron(n_I,eye(2)) + kron(eye(2), ...
        n_I)),[2 2 2 2]);

    for i =(1:length(ff_sp))
        % gg_sp(i) = 0, no on-site bath term in any site of the chain
        Hc{i} = ff_sp(i)'*contract(F,3,2,permute(conj(F),[3 2 1]),3,2);
        Hc{i} = Hc{i} + ff_sp(i)*contract(permute(conj(F), ...
            [3 2 1]),3,2,F,3,2);
    end
    % Now we take into account U and epsd for the impurity.

    % First the site that acts on the first spin up bath and spin up imp
    Hc{N} = Hc{N} + permute(reshape(epsd*kron(eye(2), ...
        n_I),[2 2 2 2]),[1 3 2 4]);

    % Now actiong on both impurities (spin up-down)
    Hc{N+1} = Hc{N+1} + permute(reshape(U*kron(n_I,n_I) ...
        ,[2 2 2 2]),[1 3 2 4]);
    Hc{N+1} = Hc{N+1} + permute(Tensor_N1,[1 3 2 4]);

    % Now acting on the spin-up imp and first spin-up bath site
    Hc{N+2} = Hc{N+2} + permute(reshape(epsd*kron(n_I,eye(2)), ...
        [2 2 2 2]),[1 3 2 4]);
end

% Construct the Hamiltonian as two site operators for Star geometry (U=0)
Hs = cell(Nsite-1,1);

if U == 0
    for i =(1:length(ff))
        Tensor_E1 = reshape(estar(i)*kron(eye(2),n_I),[2 2 2 2]);
        Tensor_E1 = permute(Tensor_E1,[1 3 2 4]);
        Hs{i} = vstar(i)*contract(F,3,2,permute(conj(F),[3 2 1]),3,2);
        Hs{i} = Hs{i} + vstar(i)*contract(permute(conj(F), ...
            [3 2 1]),3,2,F,3,2);
        Hs{i} = Hs{i} + Tensor_E1;
    end
end


%% Call time-evolution functions

if U == 0
   [Mf,Mr,EE,~,Green] = tDMRG_chain(M,Hc,Nkeep,dt/2,tmax/2);
   [Mfs,Mfr,EEs,~,GreenS] = tDMRG_star(M,Hs,Nkeep,dt/2,tmax/2);
else
   [Mfref,Mrref,EEref,~,Greenref] = tDMRG_chain(Mref,Hc,Nkeep,dt/4,tmax/2);
   [Mf,Mr,EE,~,Green] = tDMRG_chain(M,Hc,ceil(Nkeep/4),dt/2,tmax/2);
end



%% Green function calculation

if U == 0
    for iT = (1:length(Tstep)-1)
        Green(iT) = Green(iT)*exp(1i*dt*E0*iT); 
        GreenS(iT) = GreenS(iT)*exp(1i*dt*E0*iT);
    end
    % Including t = 0
    Green = [1;Green];
    GreenS = [1;GreenS];
else
    Greenref = Greenref(2:2:end);
    for iT = (1:length(Tstep)-1)
        Greenref(iT) = Greenref(iT)*exp(1i*dt*E0*iT);
        Green(iT) = Green(iT)*exp(1i*dt*E0*iT);        
    end
    % Including t = 0
    Greenref = [1;Greenref];
    Green = [1;Green];
end

%% Exact greater Green's function - only for U=0

if U == 0
    % Chain Geometry
    M = canonForm(M,Nsite,'-k');
    
    Tstep = 0:dt:tmax;
    [U1,E1] = eig((Hchain+Hchain')/2);
    Ud = U1';
    Gexact = zeros(length(Tstep),1);
    Ed = diag(E1);%*2;  
    Nmk = zeros(Nsite,1);  % N_k
    Ak2 = zeros(Nsite,1);  % abs(a_k)^2
    Cop = cell(Nsite,1);  % Creation MPOs.
    Aop = cell(Nsite,1);  % Annihilation MPOs.
    Fan = squeeze(F);
    Fdag = squeeze(F)';
    Z = squeeze(Z);

    % We compute the number of particles per site in the non-interacting
    % hamiltonian eigenbasis
    for itk = 1:(Nsite)

        % Creation operator MPO.
        Cop{1} = zeros(1,2,2,2);
        Cop{1}(1,:,1,:) = Ud(itk,1)*Fdag;
        Cop{1}(1,:,2,:) = eye(2);
        for itj = 2:(Nsite-1)
            Cop{itj} = zeros(2,2,2,2);
            Cop{itj}(1,:,1,:) = Z;
            Cop{itj}(2,:,1,:) = Ud(itk,itj)*Fdag;
            Cop{itj}(2,:,2,:) = eye(2);
        end
        Cop{Nsite} = zeros(2,2,1,2);
        Cop{Nsite}(1,:,1,:) = Z;
        Cop{Nsite}(2,:,1,:) = Ud(itk,Nsite)*Fdag;

        % Annhilation operator MPO.
        Aop{1} = zeros(1,2,2,2);
        Aop{1}(1,:,1,:) = U1(1,itk)*Fan;
        Aop{1}(1,:,2,:) = eye(2);
        for iti = 2:(Nsite-1)
            Aop{iti} = zeros(2,2,2,2);
            Aop{iti}(1,:,1,:) = Z;
            Aop{iti}(2,:,1,:) = U1(iti,itk)*Fan;
            Aop{iti}(2,:,2,:) = eye(2);
        end
        Aop{Nsite} = zeros(2,2,1,2);
        Aop{Nsite}(1,:,1,:) = Z;
        Aop{Nsite}(2,:,1,:) = U1(Nsite,itk)*Fan;

        % Apply the number operator to the initial state M.
        NM = contract(Aop{1},4,4,M{1},3,2);
        NM = contract(Cop{1},4,4,NM,5,2);
        NM = contract(conj(M{1}),3,2,NM,7,2);
        NM = permute(NM,[2,4,6,8,1,3,5,7]);
        for j = 2:(Nsite)
            NM = contract(NM,4,4,M{j},3,1);
            NM = contract(Aop{j},4,[1,4],NM,5,[3,4]);
            NM = contract(Cop{j},4,[1,4],NM,5,[4,1]);
            NM = contract(conj(M{j}),3,[1,2],NM,5,[4,1]);
        end
        Nmk(itk) = NM;%*2;
        %Ak2(itk) = abs(Ud(itk,1)*U1(1,itk));
        Ak2(itk) = abs(Ud(1,itk)*U1(itk,1));
    end
    NumberFermions_chain = sum(Nmk); % Total number of particles
    for itT = 1:numel(Tstep)
        Gexact(itT) = sum(Ak2.*exp(-1i*((Ed.*Nmk))*Tstep(itT)));
    end

    % Star Geometry
    
    Tstep = 0:dt:tmax;
    [U1,E1] = eig((Hstar+Hstar')/2);
    Ud = U1';
    GexactS = zeros(length(Tstep),1);
    Eds = diag(E1);%*2;  
    Nmk = zeros(Nsite,1);  % N_k
    Ak2 = zeros(Nsite,1);  % abs(a_k)^2
    Cop = cell(Nsite,1);  % Creation MPOs.
    Aop = cell(Nsite,1);  % Annihilation MPOs.
    
    for itk = 1:(Nsite)

        % Creation operator MPO.
        Cop{1} = zeros(1,2,2,2);
        Cop{1}(1,:,1,:) = Ud(itk,1)*Fdag;
        Cop{1}(1,:,2,:) = eye(2);
        for itj = 2:(Nsite-1)
            Cop{itj} = zeros(2,2,2,2);
            Cop{itj}(1,:,1,:) = Z;
            Cop{itj}(2,:,1,:) = Ud(itk,itj)*Fdag;
            Cop{itj}(2,:,2,:) = eye(2);
        end
        Cop{Nsite} = zeros(2,2,1,2);
        Cop{Nsite}(1,:,1,:) = Z;
        Cop{Nsite}(2,:,1,:) = Ud(itk,Nsite)*Fdag;

        % Annhilation operator MPO.
        Aop{1} = zeros(1,2,2,2);
        Aop{1}(1,:,1,:) = U1(1,itk)*Fan;
        Aop{1}(1,:,2,:) = eye(2);
        for iti = 2:(Nsite-1)
            Aop{iti} = zeros(2,2,2,2);
            Aop{iti}(1,:,1,:) = Z;
            Aop{iti}(2,:,1,:) = U1(iti,itk)*Fan;
            Aop{iti}(2,:,2,:) = eye(2);
        end
        Aop{Nsite} = zeros(2,2,1,2);
        Aop{Nsite}(1,:,1,:) = Z;
        Aop{Nsite}(2,:,1,:) = U1(Nsite,itk)*Fan;

        % Apply the number operator to the initial state M.
        NM = contract(Aop{1},4,4,M{1},3,2);
        NM = contract(Cop{1},4,4,NM,5,2);
        NM = contract(conj(M{1}),3,2,NM,7,2);
        NM = permute(NM,[2,4,6,8,1,3,5,7]);
        for j = 2:(Nsite)
            NM = contract(NM,4,4,M{j},3,1);
            NM = contract(Aop{j},4,[1,4],NM,5,[3,4]);
            NM = contract(Cop{j},4,[1,4],NM,5,[4,1]);
            NM = contract(conj(M{j}),3,[1,2],NM,5,[4,1]);
        end
        Nmk(itk) = NM;%*2;
        %Ak2(itk) = abs(Ud(itk,1)*U1(1,itk));
        Ak2(itk) = abs(Ud(1,itk)*U1(itk,1));
    end
    NumberFermions_star = sum(Nmk); % Total number of particles
    for itT = 1:numel(Tstep)
        GexactS(itT) = sum(Ak2.*exp(-1i*((Eds.*Nmk-E0))*Tstep(itT)));
    end

    
end
%% Plots

choice = menu('Do you want the Time Evolution plots?','No','Yes');
if choice==2 || choice==0
    if U == 0  % Non-Interacting Case (U=0)
        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 375 410]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [975 310]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 375 310]);
        set(gca,'units','points');
        imagesc((1:Nsite),Tstep,EE);
        xlabel('site $\ell$','Interpreter','latex','fontsize',18);
        ylabel('time $t$','Interpreter','latex','fontsize',18);
        colormap(turbo)
        colorbar('eastoutside');
        set(gca,'innerposition',[70,60,220,310]);
        set(gca,'FontSize',18);
        title('$ (S_\mathrm{bond}(t) )_\mathrm{chain}^{(U=0)}$', ...
            'interpreter','latex');

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 375 410]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [975 310]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 375 310]);
        set(gca,'units','points');
        imagesc((1:Nsite),Tstep,EEs);
        xlabel('site $\ell$','Interpreter','latex','fontsize',18);
        ylabel('time $t$','Interpreter','latex','fontsize',18);
        colormap(turbo)
        colorbar('eastoutside');
        set(gca,'innerposition',[70,60,220,310]);
        set(gca,'FontSize',18);
        title('$ (S_\mathrm{bond}(t) )_\mathrm{star}^{(U=0)}$', ...
            'interpreter','latex');

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        semilogy(Tstep,abs(real(Green)-real(Gexact)),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$| \mathcal{R} G^>_\mathrm{DMRG}(t) - ' ...
            ' \mathcal{R} G^>_\mathrm{exact}(t)|_\mathrm{Chain}, (U=0)$'], ...
            'Interpreter','latex','fontsize',18,'location','south');
        grid on
        
        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        plot(Tstep,real(Green),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$\mathcal{R} G^>_\mathrm{DMRG}(t)' ...
            ', (U=0),$ Chain Geometry'], ...
            'Interpreter','latex','fontsize',18,'location','north');
        grid on
        
        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        plot(Tstep,real(Gexact),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$\mathcal{R} G^>_\mathrm{exact}(t)' ...
            ', (U=0),$ Chain Geometry'], ...
            'Interpreter','latex','fontsize',18,'location','north');
        grid on

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        semilogy(Tstep,abs(real(GreenS)-real(GexactS)),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$| \mathcal{R} G^>_\mathrm{DMRG}(t) - ' ...
            ' \mathcal{R} G^>_\mathrm{exact}(t)|_\mathrm{Star}, (U=0)$'], ...
            'Interpreter','latex','fontsize',18,'location','south');
        grid on

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        plot(Tstep,real(GreenS),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$\mathcal{R} G^>_\mathrm{DMRG}(t)' ...
            ', (U=0),$ Star Geometry'], ...
            'Interpreter','latex','fontsize',18,'location','north');
        grid on

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        plot(Tstep,real(GexactS),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$\mathcal{R} G^>_\mathrm{exact}(t)' ...
            ', (U=0),$ Star Geometry'], ...
            'Interpreter','latex','fontsize',18,'location','north');
        grid on
      
    else       % Interacting Case (U=1)
        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        plot(Tstep,real(Greenref),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$\mathcal{R} G^>_\mathrm{DMRG}(t)' ...
            ', (U=1),$ Chain Geometry'], ...
            'Interpreter','latex','fontsize',18,'location','north');
        grid on

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 375 410]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [975 310]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 375 310]);
        set(gca,'units','points');
        imagesc((1:Nsite),Tstep,EEref);
        xlabel('site $\ell$','Interpreter','latex','fontsize',18);
        ylabel('time $t$','Interpreter','latex','fontsize',18);
        colormap(turbo)
        colorbar('eastoutside');
        set(gca,'innerposition',[70,50,220,310]);
        set(gca,'FontSize',18);
        title('$ (S_\mathrm{bond}(t) )_\mathrm{chain}^{(U=1)}$', ...
            'interpreter','latex');

        figure
        set(gcf,'units','points');
        set(gcf,'position',[0 0 700 450]);
        set(gcf,'paperunits', 'points');
        set(gcf,'papersize', [700 450]);
        set(gcf,'paperpositionmode', 'manual');
        set(gcf,'paperposition', [0 0 700 450]);
        set(gca,'units','points');
        set(gca,'innerposition',[120,50,500,350]);
        semilogy(Tstep,abs(real(Green)-real(Greenref)),'b-')
        xlabel('$t$','Interpreter','latex','fontsize',18);
        set(gca,'FontSize',18);
        set(gca,'yminortick','on');
        set(gca,'ticklength',[0.04 0.1]);
        legend(['$| \mathcal{R} G^>_\mathrm{DMRG}(t) - ' ...
            ' \mathcal{R} G^>_\mathrm{DMRG_{ref}}(t)|, (U=0)$'], ...
            'Interpreter','latex','fontsize',18,'location','south');
        grid on
    end
end


%##########################################################################
%% Functions

% Functions used in the Main code:

%   New functions:
%   - doLinDisc : linear discretization of the bath.
%   - tDMRG_chain : tDMRG algorithm for chain geometry.
%   - tDMRG_star : tDMRG algorithm for star geometry following Ref.[1].

%   Code Repository:
%   - contract
%   - DMRG_2site
%   - canonForm
%   - updateLeft
%   - svdTr
%   - getLocalSpace

%##########################################################################

%% doLinDisc.m
function [ff,gg,estarimp,vstarimp] = doLinDisc(ozin,N)
%{
< Description >
[ff,gg,estarimp,vstarimp] = doLinDisc(ozin,N)
This function performs a linear bath discretization of the hybridization
function used in Ref.[1] (semicircular), and maps the resulting 
star-geometry Hamiltonian onto the chain-geometry Hamiltonian, via Lanczos 
Tridiagonalization. The output 'ff' and 'gg'describe the hopping amplitudes 
and the on-site energies of the chain.

< Input >
ozin : [numeric vector] Linear discretization of the omega axis describing
       the hybridization function.
N : [numeric] Number of bath sites.

< Output >
ff, gg : [numeric vectors] Hopping amplitudes and on-site energies of
        the chain, respectively. The hopping amplitudes correspond
        to the superdiagonals [diag(..,+1) and diag(..,-1)] of the
        tridiagonal matrix representation of a single-particle Hamiltonian;
        the on-site energies correspond to the diagonals of the tridiagonal
        matrix.
estarimp : [numeric vector] On-site energies including including 0 at first
            site to take into account the impurity. Star Geometry
vstarimp : [numeric vector] Hopping parameter including 0 at first site to
            take into account the impurity. Star Geometry
%}

% ∆e
deltaE = ozin(2)-ozin(1);
estar = zeros(1,N);

% N equidistant intervals Ik of size ∆e. 
% estar_k = min{Ik} + ∆e/2
for i=(1:N)
    estar(i) = min(ozin(i),ozin(i+1)) + deltaE/2;
end

%{
for i=(1:N)
    a = [ozin(i) ozin(i+1)];
    b = [RhoV2in(i) RhoV2in(i+1)];
    vstar(i) = sqrt(trapz(a,b));
end
vstar
%}

% Hybridization function integral ∫√(1-w^2)dw
vstar = (ozin.*sqrt(1-ozin.^2)+asin(ozin))/2; 

% Star Hopping amplitude: (1/2π)*√(vstar(2)-vstar(1))
vstar = sqrt((vstar(2:(N+1))-vstar(1:N))/(2*pi));  


%{
estarimp = [0 estar(N:-1:(N+3)/2) estar(1:(N-1)/2) estar((N+1)/2)];
vstarimp = [0 vstar(N:-1:(N+3)/2) vstar(1:(N-1)/2) vstar((N+1)/2)];
%}

% To construct the star Hamiltonian
estarimp = [0 estar(1:N)];
vstarimp = [0 vstar(1:N)];

% Construction of the Star Hamiltonian
Hstar = diag(estarimp);
Hstar(:,1) = vstarimp;
Hstar(1,:) = vstarimp;

U = zeros(size(Hstar,1),1);
U(1,1) = 1; % Initial Krylov vector

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
choice = menu('Do you want the Bath Parameters plots?','No','Yes');
if choice==2 || choice==0
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
    title('Hopping Parameters', ...
            'interpreter','latex')
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
     title('On-site Parameters', ...
            'interpreter','latex')
end
end

%##########################################################################

%% tDMRG_chain

function [Mf,Mr,EE,dw,Green] = tDMRG_chain(M,Hs,Nkeep,dt,tmax)
%{
 < Description >
 [Mf,Mr,EE,dw,Green] = tDMRG_chain(M,Hs,Nkeep,dt,tmax)
 Real time evolution using tDMRG. The 1D chain system is 
 described by the Hamiltonian H. We employed 2nd order Trotter 
 step exp(-dt/2*Hodd) * exp(-dt*Heven) * exp(-dt/2*Hodd). After acting each
 gate, the bonds are truncated. Also this function directly computes the 
 overlap of the "forward" and "reverse" time evolutions, needed to then 
 calculate the greater Green's function.

 < Input >
 M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
       defines the chain length. 
 Hs : [cell] Hamiltonian. H{n} describes the two-site interaction between 
       site n and n+1. It should satisfy numel(M) == numel(H) + 1.
       The leg convention of H{n} are as follows:
       2      4       [legs 1 and 2 are for site n;
       |      |       legs 3 and 4 are for site n+1]
      [  H{n}  ]
       |      |
       1      3
 Nkeep : [integer] Maximum bond dimension.
 dt : [numeric] Real time step size.
 tmax : [numeric] Maximum time range.

 < Output >
 Mf(Mr) : [cell] The final MPS's after real-time evolution for the two 
          evolutions needed in the greater Green function calculation.
 EE : [matrix] EE(m,n) indicates the entanglement entropy at each bond,
       after acting the full Trotter evolution. (log2 -> ebits)
 dw : [matrix] Discarded weights (i.e., the sum of the squares of the
       discarded singular values).
 Green: [vector] value of the overlap between the MPS's corresponding to 
       the "forward" and "reverse" time evolution.
%}   

tobj = tic2;

%
% Sanity Check
if length(M) ~= (length(Hs)+1)
    error("ERROR: it should be length(M) == (length(H)+1)\n");
end
%}
N = numel(M);
Nstep = ceil(tmax/dt);

EE = zeros(Nstep,N-1);
dw = zeros(size(EE));
Green = zeros(Nstep,1);

% show information
fprintf("\n tDMRG: Real-time evolution - Chain Gemotry\n");
fprintf(['# of sites = ',sprintf('%i',numel(Hs)), ...
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
        if mod(it1,2) == 1 % True when odd
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
% Variables ending with r refer to the reverse time evolution needed for
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
 < Description >
 [M,EE,dw] = tDMRG_1sweep(M,expH,Nkeep,isright)
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
        [M{it},S,V,dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
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
        [U,S,M{it+1},dw(it)] = svdTr(T,4,[1,2],Nkeep,[]);
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
 < Description >
 Computes the overlap between the forward and reverse time evolutions.

 < Input >
 Mf(Mr) : [cell] Input MPS's.
 isleft : [logical] If true (== 1), it means that the MPS M is in 
       left-canonical form. Otherwise, right-canonical form.
 dt : [numerical] Time-Step

 < Output >
 Green : [numerical] Green function at the desired time step.
%}

N = length(Mf);

% Sanity Check
if N ~= length(Mr)
    error('Something went wrong!')
end


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

%##########################################################################

%% tDMRG_star

function [Mf,Mr,EE,dw,Green] = tDMRG_star(M,Hs,Nkeep,dt,tmax)
%{
< Description >
 [Mf,Mr,EE,dw,Green] = tDMRG_star(M,Hs,Nkeep,dt,tmax)
 Real time evolution using tDMRG. The 1D star system is 
 described by the Hamiltonian H. We employed 2nd order Trotter 
 step exp(-dt/2*Hodd) * exp(-dt*Heven) * exp(-dt/2*Hodd). After acting each
 gate, the bonds are truncated. Also this function directly computes the 
 overlap of the "forward" and "reverse" time evolutions, needed to then 
 calculate the greater Green's function.

< Input >
 M : [cell] The initial state as the MPS. The length of M, i.e., numel(M),
       defines the chain length. 
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

< Output >
 Mf(Mr) : [cell] The final MPS's after real-time evolution for the two 
          evolutions needed in the greater Green function calculation.
 EE : [matrix] EE(m,n) indicates the entanglement entropy at each bond,
       after acting the full Trotter evolution. (log2 -> ebits)
 dw : [matrix] Discarded weights (i.e., the sum of the squares of the
       discarded singular values).
 Green: [vector] value of the overlap between the MPS's corresponding to 
       the "forward" and "reverse" time evolution.

%}   

tobj = tic2;

% Sanity Check
if length(M) ~= (length(Hs)+1)
    error("ERROR: it should be length(M) == (length(H)+1)");
end

N = numel(M);
Nstep = ceil(tmax/dt);


EE = zeros(Nstep,N-1);
dw = zeros(size(EE));
Green = zeros(Nstep,1);

% SWAP gate
SWAP = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 -1]; % matrix form in order to 
                                              % multiply it for exp(H).

%{
         ---(1)-*-(2)---
             |     |
              \   /
               \ / 
          [ SWAP gate ]
               / \
              /   \
             |     |
%}
                                        
% show information
fprintf("\n tDMRG: Real-time evolution - Star Gemotry\n");
fprintf(['# of sites = ',sprintf('%i',numel(Hs)), ...
    ', Nkeep = ',sprintf('%i',Nkeep),', dt = ',sprintf('%i',dt),'\n']);

% generate the unitray operator exp(-it*H) for each two-site pairs

% Forward time evolution gates for the exp(-iHt/2)
Fgate = cell(length(Hs)-1,1); % Forward gates (1st half-evolve then SWAP)
Rgate = cell(length(Hs)-1,1); % Reverse gates (1st SWAP then half-evolve)
Lgate = cell(1,1);            % Last gate (Full evolve only -  no SWAP's)

% Reverse time evolution gates for the exp(+iHt/2)
Fgater = cell(length(Hs)-1,1); % Forward gates (1st half-evolve then SWAP)
Rgater = cell(length(Hs)-1,1); % Reverse gates (1st SWAP then half-evolve)
Lgater = cell(1,1);            % Last gate (Full evolve only -  no SWAP's)


% Creation of the gates
for it1 = (1:length(Hs))
    if ~isempty(Hs{it1})
        sdim = [size(M{it1},2),size(M{it1+1},2)];
        Htmp = permute(Hs{it1},[1 3 2 4]);
        Htmp = reshape(Htmp,[sdim(1)*sdim(2) sdim(1)*sdim(2)]);
        if it1 == length(Hs) % Last Hs acts with full time-step

            ttmp = dt/2; % Should be dt and not dt/2, but for normalization
                         % purposes I apply it twice instead of only one
                         % time.

            % Forward time
            eH = expm(-1i*ttmp*Htmp);
            Lgate{1} = reshape(eH, size(Hs{it1}));
            % Reverse time
            eHr = expm(+1i*ttmp*Htmp);
            Lgater{1} = reshape(eHr,size(Hs{it1}));
        else
            ttmp = dt/2; % Half time step for all other steps
            % Forward time
            eH = expm(-1i*ttmp*Htmp);
            Fgate{it1} = reshape(SWAP*eH, size(Hs{it1}));
            Rgate{it1} = reshape(eH*SWAP, size(Hs{it1})); 
            % Reverse time
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

    [Mf,EE1,dw1] = S_sweep(Mf,Fgate,Rgate,Lgate,Nkeep);
    [Mr,~,~] = S_sweep(Mr,Fgater,Rgater,Lgater,Nkeep);

    EE(it1,:) = EE1;
    dw(it1,:) = dw1;
    Green(it1) = star_expVal(Mf,Mr);

    str = ['Sweep #',sprintf('%i/%i',[it1,Nstep])];
    disptime(str);
end

toc2(tobj,'-v');
end

% -------------------------------------------------------------------------

function [M,EE,dw] = S_sweep(M,Fgate,Rgate,Lgate,Nkeep)
%{
 < Description >
 [M,EE,dw] = S_sweep(M,Fgate,Rgate,Lgate,Nkeep)
 Does one full sweep of the tDMRG scheme for star geometry in Ref.[1]. As
 in tDMRG.m function, it does bond-by-bond compression.

 < Input >
 M : [cell] Input MPS.
 Fgate : [cell] contains the forward (left->right) gates consisting of 
        (SWAP*exp(H)) so that we first apply the Hamiltonian and then swap.
        Also it only evolves half timestep.
 Rgate : [cell] contains the reverse (right->left) gates consisting of 
        (exp(H)*SWAP) so that we first apply the swap and then Hamiltonian.
        Also it only evolves half timestep.
 Lgate : [cell] contains only the exp(H) for sites {(N-1),N} (Last step of
        the left-right sweep). We don't apply any swap gates and here we
        have a full timestep evolution.
 Nkeep : [numeric] Maximum bond dimension.
 
 < Output >
 M : [cell] MPS after applying exp(-it*H) and truncating bonds.
 EE : [numeric vector] Entanglement entropy at each bond.
 dw : [numeric vector] Discarded weights when truncating the bond
       dimensions.
%}


N = numel(M);
EE = zeros(1,N-1);
dw = zeros(1,N-1);
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

% I had to split the last step into applying Lgate two times for
% normalization purposes.

% Last step(-I) where we do NOT apply a SWAP gate
T = contract(M{(N-1)},3,3,M{N},3,1);
T = contract(Lgate{1},4,[3 4],T,4,[2 3],[3 1 2 4]);
% SVD via svdTr
[M{(N-1)},S,V,dw(N-1)] = svdTr(T,4,[1,2],Nkeep,[]);
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
M{N} = contract(DS,2,2,V,3,1);
M{N} = M{N}/norm(M{N}(:)); % to normalize the norm of MPS

% Last step(-II) where we do NOT apply a SWAP gate
T = contract(M{(N-1)},3,3,M{N},3,1);
T = contract(Lgate{1},4,[3 4],T,4,[2 3],[3 1 2 4]);
% SVD via svdTr
[U,S,M{N},dw(N-1)] = svdTr(T,4,[1,2],Nkeep,[]);
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
% update M{end-1}
DS = diag(S);
M{N-1} = contract(U,3,3,DS,2,1);

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

function Green = star_expVal(Mf,Mr)
%{
 < Description >
 Computes the overlap between the forward and reverse time evolutions.

 < Input >
 Mf(Mr) : [cell] Input MPS's.

 < Output >
 Green : [numerical] Overlap of MPS's at a given time step.
%}

N = length(Mf);

% Sanity Check
if N ~= length(Mr)
    error('Something went wrong!')
end

%{
%We put the MPS's in right canonical form
Mf = canonForm(Mf,N,'-k');
Mf = canonForm(Mr,N,'-k');
%}

MM = []; % contraction of bra/ket tensors
for itN = (N:-1:1)
    Tf = permute(Mf{itN},[3 2 1]); % permute left<->right to use uLeft                                         
    Tr = permute(Mr{itN},[3 2 1]); % permute left<->right to use uLeft
    MM = updateLeft(MM,2,Tr,[],[],Tf);
end

Green = MM;

end


