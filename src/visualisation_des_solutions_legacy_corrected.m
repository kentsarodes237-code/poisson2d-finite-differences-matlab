%% ============================================================
%  CAS N = 4
%  Visualisation de la solution analytique et numérique
% ============================================================

clear;
clc;
close all;


%% 1. Choix de N et du pas h

N = 4;

h = 1/(N+1);

disp('Valeur du pas h :');
disp(h);


%% 2. Définition de la solution exacte et du second membre

% Solution exacte lue sur l'énoncé :
% u(x,y) = cos(pi*(x+y)) * sin(pi*(x-y))

uexact = @(x,y) cos(pi*(x+y)) .* sin(pi*(x-y));

% Comme -Delta u = f, on a :
% f(x,y) = 2*pi^2*( sin(2*pi*x) - sin(2*pi*y) )

f = @(x,y) 2*pi^2 * ( sin(2*pi*x) - sin(2*pi*y) );


%% 3. Construction des points intérieurs

% Pour N = 4 :
% x_int = [0.2 0.4 0.6 0.8]
% y_int = [0.2 0.4 0.6 0.8]

x_int = h:h:1-h;
y_int = h:h:1-h;

disp('Points interieurs en x :');
disp(x_int);

disp('Points interieurs en y :');
disp(y_int);


%% 4. Maillage intérieur avec ndgrid

% X_int(i,j) contient la coordonnée x du point intérieur (i,j)
% Y_int(i,j) contient la coordonnée y du point intérieur (i,j)

[X_int, Y_int] = ndgrid(x_int, y_int);


%% 5. Construction de la matrice A

% La matrice A correspond au schéma :
%
% 4u(i,j) - u(i-1,j) - u(i+1,j) - u(i,j-1) - u(i,j+1)
% = h^2 f(x_i,y_j)

e = ones(N,1);

% Matrice tridiagonale T
T = spdiags([-e 4*e -e], [-1 0 1], N, N);

% Matrice identité
I = speye(N);

% Matrice de couplage entre les lignes du maillage
B = spdiags([-e -e], [-1 1], N, N);

% Matrice globale A
A = kron(speye(N), T) + kron(B, I);


%% 6. Construction du second membre b

% On prépare d'abord une matrice F de taille N x N
% Elle contiendra h^2*f(x_i,y_j) plus les contributions des bords.

F = zeros(N,N);

for j = 1:N
    for i = 1:N
        
        x = X_int(i,j);
        y = Y_int(i,j);
        
        % Terme principal h^2*f(x_i,y_j)
        F(i,j) = h^2 * f(x,y);
        
        
        % Contribution du bord gauche : x = 0
        if i == 1
            F(i,j) = F(i,j) + uexact(0,y);
        end
        
        % Contribution du bord droit : x = 1
        if i == N
            F(i,j) = F(i,j) + uexact(1,y);
        end
        
        % Contribution du bord inférieur : y = 0
        if j == 1
            F(i,j) = F(i,j) + uexact(x,0);
        end
        
        % Contribution du bord supérieur : y = 1
        if j == N
            F(i,j) = F(i,j) + uexact(x,1);
        end
        
    end
end


%% 7. Transformation de F en vecteur b

% Le vecteur b doit être rangé dans le même ordre que les inconnues U.
%
% On choisit l'ordre :
% U = [u(1,1), u(2,1), u(3,1), u(4,1),
%      u(1,2), u(2,2), ..., u(4,4)]^T

b = zeros(N*N,1);

k = 1;

for j = 1:N
    for i = 1:N
        
        b(k) = F(i,j);
        
        k = k + 1;
        
    end
end


%% 8. Résolution du système linéaire

U = A\b;

disp('Vecteur solution numerique U :');
disp(U);


%% 9. Construction du maillage complet avec les bords

% Pour N = 4 :
% x_full = [0 0.2 0.4 0.6 0.8 1]
% y_full = [0 0.2 0.4 0.6 0.8 1]

x_full = 0:h:1;
y_full = 0:h:1;

[X_full, Y_full] = ndgrid(x_full, y_full);


%% 10. Calcul de la solution exacte sur tout le domaine

U_exact_full = zeros(N+2, N+2);

for j = 1:N+2
    for i = 1:N+2
        
        x = X_full(i,j);
        y = Y_full(i,j);
        
        U_exact_full(i,j) = uexact(x,y);
        
    end
end


%% 11. Construction de la solution numérique complète

% On commence par mettre la solution exacte partout.
% Cela permet de remplir correctement les bords.
% Ensuite, on remplacera seulement les points intérieurs par la solution numérique.

U_num_full = U_exact_full;


%% 12. Placement des valeurs numériques intérieures

k = 1;

for j = 1:N
    for i = 1:N
        
        % Dans la matrice complète :
        % l'indice 1 correspond au bord x=0 ou y=0.
        % Donc le point intérieur i devient i+1.
        
        ligne = i + 1;
        colonne = j + 1;
        
        U_num_full(ligne, colonne) = U(k);
        
        k = k + 1;
        
    end
end


%% 13. Affichage de la solution analytique

figure;

surf(X_full, Y_full, U_exact_full);

xlabel('x');
ylabel('y');
zlabel('u_{ex}(x,y)');

title('Solution analytique pour N = 4');

colorbar;
grid on;


%% 14. Affichage de la solution numérique

figure;

surf(X_full, Y_full, U_num_full);

xlabel('x');
ylabel('y');
zlabel('u_h(x,y)');

title('Solution numérique pour N = 4');

colorbar;
grid on;


%% 15. Comparaison côte à côte

figure;

subplot(1,2,1);

surf(X_full, Y_full, U_exact_full);

xlabel('x');
ylabel('y');
zlabel('u_{ex}');
title('Solution analytique');
colorbar;
grid on;


subplot(1,2,2);

surf(X_full, Y_full, U_num_full);

xlabel('x');
ylabel('y');
zlabel('u_h');
title('Solution numérique');
colorbar;
grid on;