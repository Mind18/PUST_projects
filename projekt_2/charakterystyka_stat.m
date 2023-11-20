clear all;
close all;

n=300;
U(1:n)=0;
Z(1:n)=0;
Y(1:7)=0;
Y(8:n)=1;
Ystat3D=zeros(101);
Ustat=zeros(1, 101);
Zstat=zeros(1, 101);
    

%Wyznaczanie charakterystki statycznej y(u,z)
for i=1:101
    dU=(i-1)*0.01;
    U(10:n)=dU;
    for j=1:101
        dZ=(j-1)*0.01;
        Z(10:n)=dZ;
        for k=8:n
            Y(k)=symulacja_obiektu8y_p2(U(k-6),U(k-7),Z(k-1),Z(k-2),Y(k-1),Y(k-2));
        end
        Zstat(j)=Z(n);
        Ystat3D(i,j)=Y(n);
    end
    Ustat(i)=U(n);
end

figure;
surf(Ustat,Zstat,Ystat3D); % Utworzenie wykresu 3D
xlabel('U');
ylabel('Z');
zlabel('Y');
title('Charakterystyka statyczna y(u,z)');
axes.SortMethod='ChildOrder';
export_fig('./pliki_wynikowe/zad2_y(u,z)_char_stat.pdf');


