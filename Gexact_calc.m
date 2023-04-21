dt = 0.05;
xt = 0:dt:20;
Gexact = zeros(numel(xt),1);
T = diag(gg) + diag(ff,1) + diag(ff,-1);  % Hopping matrix.
[V1,D1] = eig(T);  % Diagonalization.
Ed = diag(D1);  % Energies.
Cn2 = zeros(2*N+2,1);  % |C_n|^2
Crx = cell(2*N+2,1);  % Creation MPOs.
Anx = cell(2*N+2,1);  % Annihilation MPOs.
Vd = V1';
for itk = 1:(2*N+2)
    % First: build annihilation operator MPO.
    Anx{1} = zeros(1,2,2,2);
    Anx{1}(1,:,1,:) = Vd(itk,1)*sm;
    Anx{1}(1,:,2,:) = sz;
    for itj = 2:(2*N+1)
        Anx{itj} = zeros(2,2,2,2);
        Anx{itj}(1,:,1,:) = id;
        Anx{itj}(2,:,1,:) = Vd(itk,itj)*sm;
        Anx{itj}(2,:,2,:) = sz;
    end
    Anx{2*N+2} = zeros(2,2,1,2);
    Anx{2*N+2}(1,:,1,:) = id;
    Anx{2*N+2}(2,:,1,:) = Vd(itk,2*N+2)*sm;
    % Second: build creation operator MPO.
    Crx{1} = zeros(1,2,2,2);
    Crx{1}(1,:,1,:) = V1(1,itk)*sp;
    Crx{1}(1,:,2,:) = sz;
    for iti = 2:(2*N+1)
        Crx{iti} = zeros(2,2,2,2);
        Crx{iti}(1,:,1,:) = id;
        Crx{iti}(2,:,1,:) = V1(iti,itk)*sp;
        Crx{iti}(2,:,2,:) = sz;
    end
    Crx{2*N+2} = zeros(2,2,1,2);
    Crx{2*N+2}(1,:,1,:) = id;
    Crx{2*N+2}(2,:,1,:) = V1(2*N+2,itk)*sp;
    % Apply the number operator to Psi0.
    % % contract the left half.
    cl = contract(Anx{1},4,4,Psi0{1},3,2);
    cl = contract(Crx{1},4,4,cl,5,2);
    cl = contract(conj(Psi0{1}),3,2,cl,7,2);
    cl = permute(cl,[2,4,6,8,1,3,5,7]);
    for j = 2:(N+1)
        cl = contract(cl,4,4,Psi0{j},3,1);
        cl = contract(Anx{j},4,[1,4],cl,5,[3,4]);
        cl = contract(Crx{j},4,[1,4],cl,5,[4,1]);
        cl = contract(conj(Psi0{j}),3,[1,2],cl,5,[4,1]);
    end
    if numel(size(cl)) < 3
        cl = reshape(cl,[1,size(a)]);
    end
    % % contract the right half.
    cr = contract(Anx{2*N+2},4,4,Psi0{2*N+2},2,2);
    cr = contract(Crx{2*N+2},4,4,cr,4,2);
    cr = contract(conj(Psi0{2*N+2}),2,2,cr,6,2);
    cr = permute(cr,[1,2,4,6,3,5]);
    for j = (2*N+1):-1:(N+2)
        cr = contract(cr,4,4,Psi0{j},3,3);
        cr = contract(Anx{j},4,[3,4],cr,5,[3,5]);
        cr = contract(Crx{j},4,[3,4],cr,5,[4,2]);
        cr = contract(conj(Psi0{j}),3,[2,3],cr,5,[2,4]);
    end
    % % combine left and right halves.
    Cn2(itk) = contract(cl,4,[1,2,3,4],cr,4,[1,2,3,4]);
end
for itT = 1:numel(xt)
    Gexact(itT) = sum(Cn2.*exp(-1i*(Ed-E0)*xt(itT)));
end