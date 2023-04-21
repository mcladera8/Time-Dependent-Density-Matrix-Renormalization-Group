% Let us construct the MPO, notice how the Hchain has the form of 
% SUM_{nn}(t_{ij}*C^*_{i}*C_{j} + h.c) HOWEVER we don't have to take 
% into account the spin degree of freedom due to the particle-hole symmetry
% so the Hilbert space is two dimensional (empty or occupied).
% Notice how we number of bath sites --> 2N = 2(N+1) counting the impurity

U = 0;
epsd = -U/2; 
[F,Z,I] = getLocalSpace('Fermion');

% MPO for Hff: START %
% As stated in the paper by Daniel Bauernfeind et.al., we'll be using a
% local Hilbert of dimension 2 (empty |0> or occupied |1>)
 Nsite = 2*(N+1);
 W = cell(Nsite,1); % MPO W
 % First site
 W{1} = zeros(1,2,4,2); % ordering: left bottom right top
 W{1}(1,:,2,:) = -ff(N)*squeeze(F)'; % t_l*creation
 W{1}(1,:,3,:) = -ff(N)'*squeeze(F); % t_l*annihilation
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
        W{itL2}(4,:,2,:) = ff((N+1)-itL2)*squeeze(F)'; % t_l*creation
        W{itL2}(4,:,3,:) = ff((N+1)-itL2)'*squeeze(F); % t_l*annihilation
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
 W{N+2}(2,:,1,:) = squeeze(F)'*squeeze(F); % annihilation
 W{N+2}(4,:,1,:) = epsd*squeeze(F)'*squeeze(F); % epsd*number
 W{N+2}(4,:,2,:) = ff(1)*squeeze(F)'; % t_I*creation
 W{N+2}(4,:,3,:) = ff(1)*squeeze(F); % t_I*annihilation
 W{N+2}(4,:,4,:) = I;
 for itL2 = ((N+3):(Nsite-1))
        W{itL2} = zeros(4,2,4,2); % ordering: left bottom right top
        W{itL2}(1,:,1,:) = I;
        W{itL2}(2,:,1,:) = squeeze(F); % annihilation
        W{itL2}(3,:,1,:) = squeeze(F)'; % creation
        W{itL2}(4,:,2,:) = ff(itL2-(N+1))*squeeze(F)'; % t_l*creation
        W{itL2}(4,:,3,:) = ff(itL2-(N+1))'*squeeze(F); % t_l*annihilation
        W{itL2}(4,:,4,:) = I;
end
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
 %
Nkeep = 250;
Nsite = 2*(N+1);
MPS = cell(Nsite,1);
MPS{1} = rand(1,2,Nkeep);
MPS(2:Nsite-1) = {rand(Nkeep,2,Nkeep)};
MPS{Nsite} = rand(Nkeep,2,1);
 %}

%Hs,Minit,alpha,Nkeep,Nsweep
[M,E0,~] = DMRG_2site(W,MPS,3,250,100,'Econv',1e-12);

%%



% My final energy is E0 = -15.08219