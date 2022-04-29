%% =======================================================================%
% Description
% -----------
%  This program uses MPO � MPS method to caluclate the range-3 local correlation
%  matrix of the Fermi-Hubbard model.
%
%  Author : Lei Pan, panlei@mail.tsinghua.edu.cn
%  Created: 2022-04-28 19:28
%% =======================================================================%
clear all;
clc
L=256;
N=128;

M_matrix=zeros(4*(N+1),4*(N+1));
%M_matrix=cell{1,16};%zeros(4*(N+1),4*(N+1));
%m1=[-1,0;0,1];kron(eye(4),kron(m1,eye(4)))

m1=diag((-1).^(1:N+1));
m2=diag((-1).^(1:N),1);
eta_dag=zeros(4); eta_dag(4,1)=1;
eta=eta_dag';

Spin{1}=eye(4);Spin{2}=[0,1,0,0;1,0,0,0;0,0,0,0;0,0,0,0]; Spin{3}=[0,-1i,0,0;1i,0,0,0;0,0,0,0;0,0,0,0];
Spin{4}=[1,0,0,0;0,-1,0,0;0,0,0,0;0,0,0,0]; Spin{5}=[0,0,1,0;0,0,0,0;1,0,0,0;0,0,0,0];
Spin{6}=[0,0,-1i,0;0,0,0,0;1i,0,0,0;0,0,0,0]; Spin{7}=[0,0,0,0;0,0,1,0;0,1,0,0;0,0,0,0];
Spin{8}=[0,0,0,0;0,0,-1i,0;0,1i,0,0;0,0,0,0]; Spin{9}=[1,0,0,0;0,1,0,0;0,0,-2,0;0,0,0,0]/sqrt(3);
Spin{10}=[0,0,0,1;0,0,0,0;0,0,0,0;1,0,0,0];
Spin{11}=[0,0,0,-1i;0,0,0,0;0,0,0,0;1i,0,0,0];
Spin{12}=[0,0,0,0;0,0,0,1;0,0,0,0;0,1,0,0];
Spin{13}=[0,0,0,0;0,0,0,-1i;0,0,0,0;0,1i,0,0];
Spin{14}=[0,0,0,0;0,0,0,0;0,0,0,1;0,0,1,0];
Spin{15}=[0,0,0,0;0,0,0,0;0,0,0,-1i;0,0,1i,0];
Spin{16}=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,-3]/sqrt(6);

index1=1:4:4*N+1;index2=2:4:4*N+2;index3=3:4:4*N+3;index4=4:4:4*N+4;
%M_L=kron(m1,eye(4))+kron(m2,eta);
%M_matrixLR=cell{16,16};
mm=0;
for ii=1:16
%     for jj=1:16
%         for kk=1:16
   M_matrix1{ii}=sparse(kron(m1,Spin{ii}*eye(4))+kron(m2,Spin{ii}*eta_dag));
   M_matrix2{ii}=sparse(kron(m1,eye(4)*Spin{ii})+kron(m2,eta*Spin{ii}));

end

for jj=1:16
    for kk=1:16

%     M_jj=kron(m1,Spin{jj}*eye(4))+kron(m2,Spin{jj}*eta_dag);
%     M_kk=kron(m1,Spin{kk}*eye(4))+kron(m2,Spin{kk}*eta_dag);
%     M=M_R{ii};

    M_jj=M_matrix1{jj};
    M_kk=M_matrix2{kk};  %����û�и����������M_kk=conj(M_matrix2{kk});
    M_matrixLR{kk,jj}=sparse(kron(M_kk(index1,index1),M_jj(index1,index1)))+sparse(kron(M_kk(index1,index2),M_jj(index2,index1)))...
              +sparse(kron(M_kk(index1,index3),M_jj(index3,index1)))+sparse(kron(M_kk(index1,index4),M_jj(index4,index1)));

%     M_matrixL{ii}=sparse(kron(M_jj(index1,index1),M_jj(index1,index1)))+sparse(kron(M_ii(index2,index1),M_ii(index2,index1)))...
%               +sparse(kron(M_jj(index3,index1),M_ii(index3,index1)))+kron(M_ii(index4,index1),M_ii(index4,index1));
    end
end

M1=sparse(kron(m1,eye(4)))+sparse(kron(m2,eta_dag));
M1_dag=sparse(kron(m1,eye(4)))+sparse(kron(m2,eta));


b_M1=M1(1:4,:); M1_b=M1(:,end-3:end);

b_M1=zeros(1,N+1); b_M1(1,1)=1;
M1_b=zeros(N+1,1); M1_b(end,1)=1;

% M_matrix=sparse(kron(M1(index1,index1),M1(index1,index1)))+sparse(kron(M1(index2,index1),M1(index2,index1)))...
%               +sparse(kron(M1(index3,index1),M1(index3,index1)))+sparse(kron(M1(index4,index1),M1(index4,index1)));

M_matrix=sparse(kron(M1_dag(index1,index1),M1(index1,index1)))+sparse(kron(M1_dag(index1,index2),M1(index2,index1)))...
              +sparse(kron(M1_dag(index1,index3),M1(index3,index1)))+sparse(kron(M1_dag(index1,index4),M1(index4,index1)));

% Site1b=sparse(kron(b_M1(1,index1),b_M1(1,index1)))+sparse(kron(b_M1(2,index1),b_M1(2,index1)))...
%             +sparse(kron(b_M1(3,index1),b_M1(3,index1)))+sparse(kron(b_M1(4,index1),b_M1(4,index1)));

Site1b=kron(b_M1,b_M1);

% SiteLb=sparse(kron(M1_b(index1,1),M1_b(index1,1)))+sparse(kron(M1_b(index1,2),M1_b(index1,2)))...
%             +sparse(kron(M1_b(index1,3),M1_b(index1,3)))+sparse(kron(M1_b(index1,4),M1_b(index1,4)));

SiteLb=kron(M1_b,M1_b);

        Site1b*(M_matrix)^(L)*SiteLb/nchoosek(L,N)
        Site1b*(M_matrix)^(L-1)*SiteLb
        Site1b*(M_matrixLR{1,1})^(L-2)*SiteLb
        Site1b*SiteLb
   M11=M_matrixLR{1,1};
   Xishu=nchoosek(L,N);


   ML_3=M11^(L-3);
Corr=zeros(16^3-1,16^3-1); Li=cell(1,4095);
for jj=1:16^3-1
    aa=floor(jj/256)+1;
        bb=floor((jj-(aa-1)*256)/16)+1;
        cc=jj-(aa-1)*256-(bb-1)*16+1;
   Li{jj}=sparse(Site1b*M_matrixLR{1,aa}*M_matrixLR{1,bb}*M_matrixLR{1,cc}*ML_3*SiteLb)/Xishu;

end

for mm=1:16^3-1
    mm
     tic

     aa=floor(mm/256)+1;
        bb=floor((mm-(aa-1)*256)/16)+1;
        cc=mm-(aa-1)*256-(bb-1)*16+1;

    for nn=mm:16^3-1

        nn;
        dd=floor(nn/256)+1;
        ee=floor((nn-(dd-1)*256)/16)+1;
        ff=nn-(dd-1)*256-(ee-1)*16+1;


M1aa=M_matrixLR{1,aa}; M1bb=M_matrixLR{1,bb}; M1cc=M_matrixLR{1,cc};
M1dd=M_matrixLR{1,dd}; M1ee=M_matrixLR{1,ee}; M1ff=M_matrixLR{1,ff};
M1ad=M_matrixLR{aa,dd};  M1be=M_matrixLR{bb,ee}; M1cf=M_matrixLR{cc,ff};

M1da=M_matrixLR{dd,aa};  M1eb=M_matrixLR{ee,bb}; M1fc=M_matrixLR{ff,cc};

% Li=Site1b*M_matrixLR{1,aa}*M_matrixLR{1,bb}*M_matrixLR{1,cc}*M_matrixLR{1,1}^(L-3)*SiteLb;
% Lj=Site1b*M_matrixLR{1,dd}*M_matrixLR{1,ee}*M_matrixLR{1,ff}*M_matrixLR{1,1}^(L-3)*SiteLb;
% Corr(mm+1,nn+1)=Site1b*M_matrixLR{aa,dd}*M_matrixLR{bb,ee}*M_matrixLR{cc,ff}*M_matrixLR{1,1}^(L-3)*SiteLb-Li*Lj;


% M_matrixLR{1,aa}^L*M_matrixLR{1,bb}^L;



%Lj=Site1b*M1dd*M1ee*M1ff*ML_3*SiteLb/Xishu;
Corr(mm,nn)=1/2*(Site1b*M1ad*M1be*M1cf*ML_3*SiteLb+Site1b*M1da*M1eb*M1fc*ML_3*SiteLb)/Xishu-Li{mm}*Li{nn};
Corr(nn,mm)=Corr(mm,nn);
%Corr(mm,nn)

    end
    toc
end


[Vec,EE]=eig(Corr);
[EEE,JJJ]=sort(diag(EE),'ascend');

EEEE=abs(EEE);

N_zero=length(find(EEE<10^(-13)));
