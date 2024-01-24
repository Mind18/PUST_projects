import symulacja_obiektu8y_p4.*
% clear all; close all;
% load('s.mat');
wykres=1;
D=400;
nu=4;
ny=3;

s11 = s{1, 1};
s12 = s{1, 2};
s13 = s{1, 3};
s14 = s{1, 4};
s21 = s{2, 1};
s22 = s{2, 2};
s23 = s{2, 3};
s24 = s{2, 4};
s31 = s{3, 1};
s32 = s{3, 2};
s33 = s{3, 3};
s34 = s{3, 4};

%% Macierz odpowiedzi skokowych
S{900}={};
for i=1:900
    S(i)={[s11(i) s12(i) s13(i) s14(i);...
           s21(i) s22(i) s23(i) s24(i);...
           s31(i) s32(i) s33(i) s34(i)]};
end
    
%% Parametry dobrane eksperymentalnie
N=D; Nu=D;
lambda1=1; lambda2=1; lambda3=1; lambda4=1;
psi1=1; psi2=1; psi3=1;

%% Punkty pracy
U1pp=0;U2pp=0;U3pp=0;U4pp=0;
Y1pp=0;Y2pp=0;Y3pp=0;
%% Wyliczenia
Mp{N, D-1}={};
M{N, Nu}={};
size(Nu)=0;
size2(N)=0;
L{Nu, Nu}={};
Psi{N, N}={};
A{Nu, Nu}={};
%Macierz Mp
for i=1:(D-1)
    for j=1:N
        Mp{j,i}=S{j+i}-S{i};
    end
end

%Macierz M i Lambd
for i=1:Nu
    M(i:N,i)=S(1:N-i+1);
    M(1:i-1,i)={[0 0 0 0; 0 0 0 0; 0 0 0 0]};
    size(i)=nu;
    for j=1:Nu
        if i==j
            L{i,j}=[lambda1 0 0 0;...
                    0 lambda2 0 0;...
                    0 0 lambda3 0;...
                    0 0 0 lambda4];
        else
            L{i,j}=[0 0 0 0;...
                    0 0 0 0;...
                    0 0 0 0;...
                    0 0 0 0;];
        end
    end
end

%Maciersz Psi
for i=1:N
    size2(i)=ny;
    for j=1:N
        if i==j
            Psi{i,j}=[psi1 0 0;...
                      0 psi2 0;...
                      0 0 psi3];
        else
            Psi{i,j}=[0 0 0;...
                      0 0 0;...
                      0 0 0];
        end
    end
end

%Macierz K
M_tmp=cell2mat(M);
Mt_tmp=M_tmp';
Psi_tmp=cell2mat(Psi);
Mt_tmp=Mt_tmp*Psi_tmp;
Temp_M = mat2cell(Mt_tmp*M_tmp,size,size);
for i=1:Nu
    for j=1:Nu
        A{i,j}=Temp_M{i,j}+L{i,j};
    end
end
A_tmp = cell2mat(A);
A_odwr = mat2cell(A_tmp^(-1),size,size);
A_odwr_tmp = cell2mat(A_odwr);

K = mat2cell(A_odwr_tmp*Mt_tmp,size,size2);

Mp_tmp = cell2mat(Mp);
K1=cell2mat(K(1,:));
Ku=K1*Mp_tmp;
Ke(1,1)=sum(K1(1,1:3:N*ny));Ke(1,2)=sum(K1(1,2:3:N*ny));Ke(1,3)=sum(K1(1,3:3:N*ny));
Ke(2,1)=sum(K1(2,1:3:N*ny));Ke(2,2)=sum(K1(2,2:3:N*ny));Ke(2,3)=sum(K1(2,3:3:N*ny));
Ke(3,1)=sum(K1(3,1:3:N*ny));Ke(3,2)=sum(K1(3,2:3:N*ny));Ke(3,3)=sum(K1(3,3:3:N*ny));
Ke(4,1)=sum(K1(4,1:3:N*ny));Ke(4,2)=sum(K1(4,2:3:N*ny));Ke(4,3)=sum(K1(4,3:3:N*ny));


dup=zeros(nu,1);
dUp_m{D-1, 1}={};
Y_zad_m{N, 1}={};
Y_m{N, 1}={};
Y0{N,1}={};
dU_m{Nu}={};
du=zeros(nu,1);
Y_zad_dmc=zeros(ny,1);
Y_dmc=zeros(ny,1);

for i=1:D-1
    dUp_m(i,1)={dup};
end
for i=1:N
    Y_zad_m(i,1)={Y_zad_dmc};
    Y_m(i,1)={Y_dmc};
    Y0(i,1)={[0;0]};
end
for i=1:Nu
    dU_m(i,1)={du};
end

n=400;
U1(1:n)=U1pp;U2(1:n)=U2pp;U3(1:n)=U3pp;U4(1:n)=U4pp;
Y1(1:n)=Y1pp;Y2(1:n)=Y2pp;Y3(1:n)=Y3pp;

Y1_zad(1:5)=Y1pp;Y1_zad(6:100)=1;Y1_zad(101:200)=1;Y1_zad(201:300)=1.5;Y1_zad(301:400)=0.5;
Y2_zad(1:5)=Y2pp;Y2_zad(6:100)=1;Y2_zad(101:200)=0.5;Y2_zad(201:300)=1;Y2_zad(301:400)=1.5;
Y3_zad(1:5)=Y3pp;Y3_zad(6:100)=1;Y3_zad(101:200)=0.5;Y3_zad(201:300)=1.5;Y3_zad(301:400)=1;

u1=U1-U1pp;u2=U2-U2pp;u3=U3-U3pp;u4=U4-U4pp;
y1_zad=Y1_zad-Y1pp;y2_zad=Y2_zad-Y2pp;y3_zad=Y3_zad-Y3pp;
y1(1:n)=0;y2(1:n)=0;y3(1:n)=0;
Ey=zeros(3,n);%Eu=zeros(1,n);
%% Symulacja
for k=5:n
    [Y1(k),Y2(k),Y3(k)]=symulacja_obiektu8y_p4(U1(k-1), U1(k-2), U1(k-3), U1(k-4),...
    U2(k-1), U2(k-2), U2(k-3), U2(k-4), U3(k-1), U3(k-2), U3(k-3), U3(k-4),...
    U4(k-1), U4(k-2), U4(k-3), U4(k-4), Y1(k-1), Y1(k-2), Y1(k-3), Y1(k-4),...
    Y2(k-1), Y2(k-2), Y2(k-3), Y2(k-4), Y3(k-1), Y3(k-2), Y3(k-3), Y3(k-4));
    y1(k)=Y1(k)-Y1pp;y2(k)=Y2(k)-Y2pp;y3(k)=Y3(k)-Y3pp;
    % Uchyb
    Ey(1,k)=y1_zad(k)-y1(k);Ey(2,k)=y2_zad(k)-y2(k);Ey(3,k)=y3_zad(k)-y3(k);
    
    K1_tmp=Ke*Ey(:,k);

    Y_dmc(1)=y1(k);Y_dmc(2)=y2(k);Y_dmc(3)=y3(k);
    for i=1:N
        Y_m(i,1)={Y_dmc};
    end
    
    Y_zad_dmc(1)=y1_zad(k);Y_zad_dmc(2)=y2_zad(k);Y_zad_dmc(3)=y3_zad(k);
    
    for i=1:N
        Y_zad_m(i,1)={Y_zad_dmc};
    end
    
    dUp_tmp = cell2mat(dUp_m);
    Ku_tmp=Ku*dUp_tmp;
    du=K1_tmp-Ku_tmp;

    for n=D-1:-1:2
      dUp_m(n)=dUp_m(n-1);
    end
   
    u1(k)=u1(k-1)+du(1);u2(k)=u2(k-1)+du(2);u3(k)=u3(k-1)+du(3);u4(k)=u4(k-1)+du(4);
    dUp_m(1,1)={du};
    
    U1(k)=u1(k)+U1pp;U2(k)=u2(k)+U2pp;U3(k)=u3(k)+U3pp;U4(k)=u4(k)+U4pp;
    
    Lambda_tmp = [lambda1 0 0 0;...
                  0 lambda2 0 0;...
                  0 0 lambda3 0;...
                  0 0 0 lambda4];
    %Eu(k) = du'*Lambda_tmp*du;
end
 
EY=norm(Ey)^2
%EU=norm(Eu)^2
%E=EY+EU
if wykres == 1
    save("Params.mat", "U1", "U2", "U3", "U4", "Y1", "Y1_zad", "Y2", "Y2_zad", "Y3", "Y3_zad");
    % rysuj_wykresy();
end

