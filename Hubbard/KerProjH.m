%% =======================================================================%
%  Description
%  -----------
%   This original Matlab implementation of the  recursive method for
%   calculating the ground-state manifold of a projector Hamiltonian.
%
%   Author: Shang Liu, sliu.phys@gmail.com
%  Created: 2022-04-28 19:31
%% =======================================================================%
function [Deg_OBC,Deg_PBC,LeftExt,RightExt]=KerProjH(P,SingleSiteDim,NSites_min,NSites_max)
    %----------------------------------------------------------------------
    % P is a positive semidefinite operator that acts on finite number of
    % sites. SingleSiteDim is the Hilbert space dimension of each single
    % site.
    % Now consider the 1D Hamiltonian H=\sum_i P_i where P_i acts on a
    % series of consecutive sites starting from the i-th site. The goal of
    % this function is to determine the degeneracy of the zero-energy
    % ground states of H, with OBC and possibly also with PBC, for syste
    %----------------------------------------------------------------------
    if nargout==1
        PBC_Spectrum=false; % No need to compute PBC spectrum in this case.
    else
        PBC_Spectrum=true;
    end

    Dim_P=size(P,1);
    Range=round(log(Dim_P)/log(SingleSiteDim));
    NSites_array=NSites_min:NSites_max;

    assert(Range>1,'We assume that the projector range is greater than 1!');

    %----------------------------------------------------------------------
    % Preallocation.
    %----------------------------------------------------------------------
    N_NSites=NSites_max-NSites_min+1;
    Deg_OBC=zeros(N_NSites,1);
    if PBC_Spectrum
        Deg_PBC=Deg_OBC;
    end

    RightExt=cell(N_NSites,1);
    LeftExt=cell(N_NSites,1);

    %----------------------------------------------------------------------
    % Compute the kernel dimension.
    %----------------------------------------------------------------------
    % OBC Hamiltonian at the minimal system size.
    Hmin_OBC=TransInvH(P,SingleSiteDim,NSites_min,'obc');

    for i_NSites=1:N_NSites

        NSites=NSites_array(i_NSites);
        %----------------------------------------------------------------------
        % Minimal system size.
        %----------------------------------------------------------------------
        if NSites==NSites_min
            H_OBC=Hmin_OBC;
            [Deg_OBC(1),GndSpace]=GetGnd(H_OBC);
            RightExt{1}=GndSpace;
            LeftExt{1}=GndSpace;

            if PBC_Spectrum
                H_PBC=Hmin_OBC;
                for j=1:(Range-1)
                    H_PBC=H_PBC+OnsiteOp(P,SingleSiteDim,NSites+1-j,NSites);
                end
                [Deg_PBC(1),~]=GetGnd(H_PBC);
            end
        %----------------------------------------------------------------------
        % Larger system sizes.
        %----------------------------------------------------------------------
        else
            %----------------------------------------------------------------------
            % OBC Left extension calculation.
            %----------------------------------------------------------------------
            if NSites<NSites_min+Range

                PreviousGndSpace=GndSpace; % Although we call the matrix 'space', it is really a particular choice of orthonormal basis.

                % For small system sizes, construct the full Hamiltonian.
                H_OBC=kron(P,speye(SingleSiteDim^(NSites-Range)));

                % Tensor the previous GndSpace with one more site.
                SubSpace=kron(speye(SingleSiteDim),PreviousGndSpace);

            else % NSites>=NSites_min+Range
                % Construct the Hamiltonian in a small subspace.
                H_OBC=kron(P,speye(Deg_OBC(i_NSites-Range)));

                SubSpace=1;
                for j=1:(Range-1)
                    SubSpace=SubSpace*LeftExt{i_NSites-Range+j};
                    SubSpace=kron(eye(SingleSiteDim),SubSpace);
                end
            end

            % Project the Hamiltonian to the smaller subspace.
            H_OBC_Sub=SubSpace'*H_OBC*SubSpace;
            H_OBC_Sub=(H_OBC_Sub+H_OBC_Sub')/2;

            [Deg_OBC(i_NSites),GndSpace_1Step]=GetGnd(H_OBC_Sub);

            % Update GndSpace.
            if NSites<NSites_min+Range % GndSpace is only needed for small system sizes.
                GndSpace=SubSpace*GndSpace_1Step; % Embed back to the full Hilbert space.
            end

            % Save the single-step embedding matrix.
            LeftExt{i_NSites}=GndSpace_1Step;

            if PBC_Spectrum
                %----------------------------------------------------------------------
                % PBC calculation.
                %----------------------------------------------------------------------
                if NSites<NSites_min+Range
                    H_PBC=sparse(SingleSiteDim^NSites,SingleSiteDim^NSites);
                    for j=1:(Range-1)
                        H_PBC=H_PBC+OnsiteOp(P,SingleSiteDim,NSites+1-j,NSites);
                    end
                    H_PBC_Sub=GndSpace'*H_PBC*GndSpace;
                    H_PBC_Sub=(H_PBC_Sub+H_PBC_Sub')/2;

                    % Be careful not to rewrite GndSpace here.
                    [Deg_PBC(i_NSites),~]=GetGnd(H_PBC_Sub);

                else
                    H_PBC_Sub=zeros(Deg_OBC(i_NSites),Deg_OBC(i_NSites));

                    for j=1:(Range-1)
                        M=kron(P,eye(Deg_OBC(i_NSites-Range)));
                        Dim1=SingleSiteDim^j;
                        Dim2=(SingleSiteDim^(Range-j))*Deg_OBC(i_NSites-Range);
                        M=SwappedOp(M,Dim1,Dim2);

                        SubSpace=RightExt{i_NSites-Range+1};
                        for iR=2:j
                            SubSpace=kron(SubSpace,eye(SingleSiteDim));
                            SubSpace=SubSpace*RightExt{i_NSites-Range+iR};
                        end

                        for iL=(j+1):Range
                            SubSpace=kron(eye(SingleSiteDim),SubSpace);
                            SubSpace=SubSpace*LeftExt{i_NSites-Range+iL};
                        end

                        H_PBC_Sub=H_PBC_Sub+SubSpace'*M*SubSpace;
                    end

                    H_PBC_Sub=(H_PBC_Sub+H_PBC_Sub')/2;

                    % Be careful not to rewrite GndSpace here.
                    [Deg_PBC(i_NSites),~]=GetGnd(H_PBC_Sub);
                end
                %----------------------------------------------------------------------
                % OBC Right extension calculation.
                %----------------------------------------------------------------------
                if NSites<NSites_min+Range

                    % PreviousGndSpace has been defined.
                    % Current GndSpace is also updated.
                    RightExt{i_NSites}=kron(PreviousGndSpace,eye(SingleSiteDim))'*GndSpace;

                else % NSites>=NSites_min+Range

                    RightExt{i_NSites}=kron(LeftExt{i_NSites-1},eye(SingleSiteDim))'*...
                        (kron(eye(SingleSiteDim),RightExt{i_NSites-1})*LeftExt{i_NSites});

                end
            end
        end
    end
end

function [Deg,GndSpace]=GetGnd(H)
    [Eigenstates,Eigenvals]=eig(full(H));
    Eigenvals=diag(Eigenvals);
    [Eigenvals,SortIndex]=sort(Eigenvals,'ascend','ComparisonMethod','real');
    Eigenstates=Eigenstates(:,SortIndex);
    GndIndices=Eigenvals<1e-12;
    Deg=sum(GndIndices);
    GndSpace=Eigenstates(:,GndIndices);
end

function Mnew=SwappedOp(M,Dim1,Dim2)
    %----------------------------------------------------------------------
    % Let A and B be two Hilbert spaces. There is an isomorphism between
    % A\otimes B and B\otimes A. Let M be an operator acting on A\otimes B
    % under the kron product representation, Mnew is the corresponding
    % operator acting on B\otimes A.
    %----------------------------------------------------------------------
    Dim=Dim1*Dim2;
    [i_row,i_col,values]=find(M);
    New_i_row=SwapIndMap(i_row,Dim1,Dim2);
    New_i_col=SwapIndMap(i_col,Dim1,Dim2);
    Mnew=sparse(New_i_row,New_i_col,values,Dim,Dim);
end


function M=OnsiteOp(M0,SingleSiteDim,SiteIndex,N_Sites)
    %----------------------------------------------------------------------
    % Put the operator M0 on a one-dimensional chain. M0 itself can act on
    % multiple number of consecutive sites. SingleSiteDim is the Hilbert
    % space dimension on each site. SiteIndex is the position of the
    % leading site in the support of M0 (after putting on the chain). We
    % label the sites from 1 to N_Sites.
    %
    % We assume M0 is a square matrix and will not check this explicitly.
    %----------------------------------------------------------------------
    if ~issparse(M0)
        M0=sparse(M0); % Convert to a sparse matrix.
    end

    Dim_M0=size(M0,1);
    StringLength=round(log(Dim_M0)/log(SingleSiteDim));

    assert((1<=SiteIndex)&&(SiteIndex<=N_Sites),'Invalid SiteIndex!');
    assert(StringLength<=N_Sites,'The operator string is too long!');

    M=kron(M0,speye(SingleSiteDim^(N_Sites-StringLength))); % Initialization.

    % Break the Hilbert space into the product of two parts.
    Dim1=SingleSiteDim^(N_Sites-SiteIndex+1);
    Dim2=SingleSiteDim^(SiteIndex-1);
    Dim=SingleSiteDim^N_Sites;

    % M=reshape(permute(reshape(M,[Dim2,Dim1,Dim2,Dim1]),[2,1,4,3]),[Dim,Dim]);
    [i_row,i_col,values]=find(M);
    New_i_row=SwapIndMap(i_row,Dim1,Dim2);
    New_i_col=SwapIndMap(i_col,Dim1,Dim2);
    M=sparse(New_i_row,New_i_col,values,Dim,Dim);
end
%{
M1=OnsiteOp(PauliMString('XZ'),2,4,4);
M2=OnsitePauliMString('XZ',4,4);
disp(max(abs(M1-M2),[],'all'));
%}


function NewIndices=SwapIndMap(Indices,m,n)
    %----------------------------------------------------------------------
    % Let A and B be two Hilbert spaces of dimensions m and n,
    % respectively. The basis for A\otimes B or
    % B\otimes A are chosen according to the canonical Kronecker product
    % convention. For example, basis vectors of A\otimes B are ordered as
    % |1>|1>,|1>|2>,...,|1>|n>,|2>|1>,|2>|2>,...,|2>|n>,...,|m>|n>.
    %
    % the input Indices is a column of basis vector indices in A\otimes B,
    % the output NewIndices is the corresponding basis vector indices in
    % B\otimes A.
    %----------------------------------------------------------------------
    [iB,iA]=ind2sub([n,m],Indices);
    NewIndices=sub2ind([m,n],iA,iB);
end
%{
P=SwapMat(2,3);
u1=kron([1.1;2.2],[pi;-3;5*1i]);
v1=kron([pi;-3;5*1i],[1.1;2.2]);
disp(v1-P*u1);
disp(v1(SwapIndMap(1:6,2,3))-u1(1:6));
%}


function H=TransInvH(M0,SingleSiteDim,N_Sites,LeadingSiteList)
    %----------------------------------------------------------------------
    % Construct a translation invariant Hamiltonian.
    % M0 is a local operator to be added, presumably repeatedly, to the
    % Hamiltonian.
    % SingleSiteDim is the Hilbert space dimension on each site.
    % N_Sites is the total number of sites.
    %
    % LeadingSiteList is a vector specifying a set of sites. For each site
    % index i in this list, there will be a term whose leading index is i.
    %
    % LeadingSiteList can also be a char string specifying the boundary
    % condition: 'obc' or 'pbc'.
    %----------------------------------------------------------------------
    if ~issparse(M0)
        M0=sparse(M0); % Convert to a sparse matrix.
    end

    Dim_M0=size(M0,1);
    TermLength=round(log(Dim_M0)/log(SingleSiteDim));
    Dim=SingleSiteDim^N_Sites; % Hilbert space dimension.

    if nargin<3
        LeadingSiteList='obc';
    end

    if ischar(LeadingSiteList)
        if strcmp(LeadingSiteList,'obc')
            LeadingSiteList=1:(N_Sites-TermLength+1);
        elseif strcmp(LeadingSiteList,'pbc')
            LeadingSiteList=1:N_Sites;
        else
            error('Invalid boundary condition!');
        end
    end

    H=sparse(Dim,Dim);
    for i_Site=LeadingSiteList
        H=H+OnsiteOp(M0,SingleSiteDim,i_Site,N_Sites);
    end
end
