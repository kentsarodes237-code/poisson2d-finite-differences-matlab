%% Poisson 2D by finite differences
% Corrected reference script for the project.
%
% Problem:
%   -Delta u = f on (0,1)x(0,1)
%   u = g on the boundary
%
% Exact solution:
%   u(x,y) = cos(pi(x+y)) sin(pi(x-y))
%
% The five-point finite-difference stencil is used and the numerical
% convergence is studied for N = 4, 8, 16, 32, 64 interior points
% per direction.

clear;
clc;
close all;

%% 1. Grid sizes
N_values = [4 8 16 32 64];
nb_tests = length(N_values);

h_values = zeros(nb_tests,1);
cond_values = zeros(nb_tests,1);
erreur_Linf_values = zeros(nb_tests,1);
erreur_L2_values = zeros(nb_tests,1);

%% 2. Exact solution and right-hand side
uexact = @(x,y) cos(pi*(x+y)) .* sin(pi*(x-y));

% Since -Delta u = f:
f = @(x,y) 2*pi^2 * (sin(2*pi*x) - sin(2*pi*y));

%% 3. Main loop
for test = 1:nb_tests

    N = N_values(test);
    h = 1/(N+1);
    h_values(test) = h;

    fprintf('\n========================================\n');
    fprintf('Calcul pour N = %d\n', N);
    fprintf('Pas h = %.10e\n', h);
    fprintf('========================================\n');

    %% Interior grid
    x_int = h:h:1-h;
    y_int = h:h:1-h;
    [X_int, Y_int] = ndgrid(x_int, y_int);

    %% Five-point finite-difference matrix
    %
    % 4u(i,j) - u(i-1,j) - u(i+1,j)
    %          - u(i,j-1) - u(i,j+1)
    % = h^2 f(x_i,y_j)

    e = ones(N,1);
    T = spdiags([-e 4*e -e], [-1 0 1], N, N);
    I = speye(N);
    B = spdiags([-e -e], [-1 1], N, N);

    % Correct sign: B already contains -1 off-diagonal entries.
    A = kron(speye(N), T) + kron(B, I);

    if N == 4
        disp('Matrice A pour N = 4 :');
        disp(full(A));
        fprintf('Verification A(1,2) = %.0f (attendu -1)\n', full(A(1,2)));
        fprintf('Verification A(1,5) = %.0f (attendu -1)\n', full(A(1,5)));
    end

    cond_values(test) = condest(A);

    %% Right-hand side including Dirichlet boundary contributions
    F = zeros(N,N);

    for j = 1:N
        for i = 1:N

            x = X_int(i,j);
            y = Y_int(i,j);

            F(i,j) = h^2 * f(x,y);

            if i == 1
                F(i,j) = F(i,j) + uexact(0,y);
            end
            if i == N
                F(i,j) = F(i,j) + uexact(1,y);
            end
            if j == 1
                F(i,j) = F(i,j) + uexact(x,0);
            end
            if j == N
                F(i,j) = F(i,j) + uexact(x,1);
            end

        end
    end

    % MATLAB column-major ordering matches the chosen unknown ordering.
    b = F(:);

    %% Linear solve
    U = A\b;
    U_num_int = reshape(U, N, N);
    U_exact_int = uexact(X_int, Y_int);

    %% Errors
    E = abs(U_num_int - U_exact_int);

    erreur_Linf_values(test) = max(E(:));

    % Discrete L2 norm scaled by the mesh area
    erreur_L2_values(test) = h * norm(E(:), 2);

    fprintf('condest(A) = %.6e\n', cond_values(test));
    fprintf('Erreur Linf = %.6e\n', erreur_Linf_values(test));
    fprintf('Erreur L2   = %.6e\n', erreur_L2_values(test));

    %% Full grid for visualization
    x_full = 0:h:1;
    y_full = 0:h:1;
    [X_full, Y_full] = ndgrid(x_full, y_full);

    U_exact_full = uexact(X_full, Y_full);
    U_num_full = U_exact_full;
    U_num_full(2:N+1, 2:N+1) = U_num_int;

    %% Exact vs numerical solution
    fig = figure('Visible','off');
    subplot(1,2,1);
    surf(X_full, Y_full, U_exact_full);
    xlabel('x'); ylabel('y'); zlabel('u_{ex}');
    title(sprintf('Solution analytique - N = %d', N));
    colorbar; grid on;

    subplot(1,2,2);
    surf(X_full, Y_full, U_num_full);
    xlabel('x'); ylabel('y'); zlabel('u_h');
    title(sprintf('Solution numerique - N = %d', N));
    colorbar; grid on;

    exportgraphics(fig, sprintf('comparaison_solutions_N%d.png',N));
    close(fig);

    %% Absolute error
    fig = figure('Visible','off');
    surf(X_int, Y_int, E);
    xlabel('x'); ylabel('y'); zlabel('|u_h-u_{ex}|');
    title(sprintf('Erreur absolue - N = %d', N));
    colorbar; grid on;
    exportgraphics(fig, sprintf('erreur_absolue_N%d.png',N));
    close(fig);

    %% Matrix sparsity
    fig = figure('Visible','off');
    spy(A);
    title(sprintf('Structure creuse de A - N = %d', N));
    exportgraphics(fig, sprintf('structure_matrice_A_N%d.png',N));
    close(fig);

end

%% 4. Summary table
resultats = table( ...
    N_values(:), ...
    h_values, ...
    cond_values, ...
    erreur_Linf_values, ...
    erreur_L2_values, ...
    'VariableNames', {'N','h','Conditionnement','Erreur_Linf','Erreur_L2'} ...
);

disp(resultats);

%% 5. Estimated convergence order
coeff_Linf = polyfit(log(h_values), log(erreur_Linf_values), 1);
coeff_L2   = polyfit(log(h_values), log(erreur_L2_values), 1);

ordre_Linf = coeff_Linf(1);
ordre_L2   = coeff_L2(1);

fprintf('\nOrdre estime en norme Linf : %.4f\n', ordre_Linf);
fprintf('Ordre estime en norme L2   : %.4f\n', ordre_L2);

fig = figure('Visible','off');
loglog(h_values, erreur_Linf_values, '-o', 'LineWidth', 1.5);
hold on;
loglog(h_values, erreur_L2_values, '-s', 'LineWidth', 1.5);
grid on;
xlabel('h');
ylabel('Erreur');
title('Evolution des erreurs en fonction du pas h');
legend( ...
    sprintf('Erreur L^\\infty, p = %.2f', ordre_Linf), ...
    sprintf('Erreur L^2, p = %.2f', ordre_L2), ...
    'Location', 'best' ...
);
exportgraphics(fig, 'courbes_loglog_erreurs.png');
close(fig);
