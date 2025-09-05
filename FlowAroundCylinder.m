close all; clear; clc

% Domain geometry
Ly = 1;
AR = 3;
Lx = AR * Ly;
Nhy = 79;
Nhx = AR * Nhy;
Ny = Nhy + 1;
Nx = Nhx + 1;
y = linspace(0, Ly, Ny);
x = linspace(0, Lx, Nx);
h = x(2) - x(1);
hq = h * h;
[X, Y] = meshgrid(x, y);

% Boundary conditions
Q = 1;
PsiN = Q * ones(1, Nx)';
PsiS = zeros(1, Nx)';
PsiW = linspace(0, Q, Ny);
PsiE = PsiW;
PsiCyl = Q / 2;

% Physical parameters
T = 30;
Re = 500;
Cou = 0.3;
beta = 0.5;
uRef = Q / Ly;
Dt = min(Cou * h / uRef, beta * Re * hq);
Nt = ceil(T / Dt) + 1;
t = linspace(0, T, Nt);
Dt = t(2) - t(1);

% Obstacle geometry
Center = [Lx/6, Ly/2];
Radius = Ly / 8;

% Domain topology and discrete Laplacian operator
k = 0;
G = zeros(Nx, Ny);
for j = 2 : Ny-1
    for i = 2 : Nx-1
        if ([x(i), y(j)] - Center) * ([x(i), y(j)] - Center)' <= Radius^2
            G(i, j) = 0;
        else
            k = k + 1;
            G(i, j) = k;
        end
    end
end
Lap = -delsq(G) ./ hq;
for i = 2 : Nx-1
    for j = 2 : Ny-1
        if G(i,j) == 0
            G(i,j) = -1;
        end
    end
end

% Domain visualization
figurePosition = [0, 50, 1640, 1640/AR];
f = figure(Position=figurePosition, Units="pixels");
D = zeros(Nx, Ny);
D(G == -1) = 1;
D(G == 0) = 1;
D = D';
spy(D)
xlabel('\it i')
ylabel('\it j')
set(gca, "YDir", "normal", "FontSize", 18)
exportgraphics(f, 'Domain.pdf', Resolution=300)

% Array initialization
zita = zeros(k, 1);
PSI = zeros(Nx, Ny);
ZITA = PSI;
U = PSI;
V = PSI;

% Boundary conditions into the arrays
U(1, :) = uRef;
U(end, :) = uRef;
PSI(1, :) = PsiW;
PSI(end, :) = PsiE;
PSI(:, 1) = PsiS;
PSI(:, end) = PsiN;
PSI(G == -1) = PsiCyl;

% Pathlines initial condition (point of injection)
xP{1} = [0, Ly/2+0.100];
xP{2} = [0, Ly/2+0.125];
xP{3} = [0, Ly/2+0.150];
xP{4} = [0, Ly/2-0.100];
xP{5} = [0, Ly/2-0.125];
xP{6} = [0, Ly/2-0.150];

% Streaklines array initialization
xS{1} = zeros(2, Nt);
xS{2} = zeros(2, Nt);
xS{3} = zeros(2, Nt);
xS{4} = zeros(2, Nt);
xS{5} = zeros(2, Nt);
xS{6} = zeros(2, Nt);
for it = 1 : Nt
    xS{1}(:, it) = xP{1};
    xS{2}(:, it) = xP{2};
    xS{3}(:, it) = xP{3};
    xS{4}(:, it) = xP{4};
    xS{5}(:, it) = xP{5};
    xS{6}(:, it) = xP{6};
end

% Graphics preparation
f = figure(Position=figurePosition, Units="pixels");
Uquiv = zeros(Nx, Ny);
Vquiv = zeros(Nx, Ny);
iq = 1 : 3 : Nx;
jq = [1 : 6 : floor(Ny/3), ...
      floor(Ny/3+3) : 3 : floor(2*Ny/3), floor(2*Ny/3+6) : 6 : Ny];
Map = [linspace(75/255, 1, 500)', linspace(90/255, 1, 500)', ...
linspace(241/255, 1, 500)'; linspace(1, 180/255, 500)', ...
linspace(1, 10/255, 500)', linspace(1, 32/255, 500)'];
DeltaPictures = round(Nt / 20, 1, "significant");

% Evolution and graphics
for it = 1 : Nt
    time = t(it);
    ZITA(G == -1) = 0;
    ZITA = ThomFormulae(ZITA, PSI, G, hq);
    ZITA = RK4(time, ZITA, ...
               @(t_, y_) ZitaRHS(t_, y_, Re, h, hq, U, V, G, Nx, Ny), ...
               Dt);
    [PSI, zita] = zita2psi(zita, ZITA, PSI, Lap, G, Nx, Ny, hq, ...
                           PsiS, PsiN, PsiW, PsiE, PsiCyl);
    [U, V] = psi2uv(PSI, U, V, G, h, Nx, Ny, uRef);
    for k = 1 : it
        xS{1}(:, k) = ...
            RK4(time, xS{1}(:, k), @(t_, y_) ...
            PathlinesRHS(t_, y_, x, y, U, V, h, h, Nx, Ny), ...
            Dt);
        xS{2}(:, k) = ...
            RK4(time, xS{2}(:, k), @(t_, y_) ...
            PathlinesRHS(t_, y_, x, y, U, V, h, h, Nx, Ny), ...
            Dt);
        xS{3}(:, k) = ...
            RK4(time, xS{3}(:, k), @(t_, y_) ...
            PathlinesRHS(t_, y_, x, y, U, V, h, h, Nx, Ny), ...
            Dt);
        xS{4}(:, k) = ...
            RK4(time, xS{4}(:, k), @(t_, y_) ...
            PathlinesRHS(t_, y_, x, y, U, V, h, h, Nx, Ny), ...
            Dt);
        xS{5}(:, k) = ...
            RK4(time, xS{5}(:, k), @(t_, y_) ...
            PathlinesRHS(t_, y_, x, y, U, V, h, h, Nx, Ny), ...
            Dt);
        xS{6}(:, k) = ...
            RK4(time, xS{6}(:, k), @(t_, y_) ...
            PathlinesRHS(t_, y_, x, y, U, V, h, h, Nx, Ny), ...
            Dt);
    end
    indices = linspace(1, it, 20000);
    if mod(it, 10) == 0
        ZITA(G == -1) = NaN;
        pcolor(x, y, ZITA');
        colormap(Map)
        shading interp
        colorbar
        clim([0.7*min(min(ZITA(G > 0))), 0.7*max(max(ZITA(G > 0)))])
        hold on
        Uquiv(iq, jq) = U(iq, jq);
        Vquiv(iq, jq) = V(iq, jq);
        quiver(X, Y, Uquiv', Vquiv', 3, Color='k')
        rectangle( ...
          Position=[Center-Radius*2.145/2, 2.145*Radius, 2.145*Radius], ...
          Curvature=[1, 1], ...
          LineWidth=2.5, FaceColor='k')
        xSplot{1} = interp1(1:it, xS{1}(1, 1:it), indices, "pchip");
        ySplot{1} = interp1(1:it, xS{1}(2, 1:it), indices, "pchip");
        xSplot{2} = interp1(1:it, xS{2}(1, 1:it), indices, "pchip");
        ySplot{2} = interp1(1:it, xS{2}(2, 1:it), indices, "pchip");
        xSplot{3} = interp1(1:it, xS{3}(1, 1:it), indices, "pchip");
        ySplot{3} = interp1(1:it, xS{3}(2, 1:it), indices, "pchip");
        xSplot{4} = interp1(1:it, xS{4}(1, 1:it), indices, "pchip");
        ySplot{4} = interp1(1:it, xS{4}(2, 1:it), indices, "pchip");
        xSplot{5} = interp1(1:it, xS{5}(1, 1:it), indices, "pchip");
        ySplot{5} = interp1(1:it, xS{5}(2, 1:it), indices, "pchip");
        xSplot{6} = interp1(1:it, xS{6}(1, 1:it), indices, "pchip");
        ySplot{6} = interp1(1:it, xS{6}(2, 1:it), indices, "pchip");
        plot(xSplot{1}, smooth(ySplot{1}), LineWidth=1.5, Color='r')
        plot(xSplot{2}, smooth(ySplot{2}), LineWidth=1.5, Color='r')
        plot(xSplot{3}, smooth(ySplot{3}), LineWidth=1.5, Color='r')
        plot(xSplot{4}, smooth(ySplot{4}), LineWidth=1.5, Color='b')
        plot(xSplot{5}, smooth(ySplot{5}), LineWidth=1.5, Color='b')
        plot(xSplot{6}, smooth(ySplot{6}), LineWidth=1.5, Color='b')
        title([ ...
            '{\it t} = ', num2str(time, 2), '. {\it i_t} = ', ...
            num2str(it), ' of ', num2str(Nt) ...
            ], FontWeight="normal")
        axis image
        axis([0, Lx, 0, Ly])
        xlabel('{\it x}')
        ylabel('{\it y}')
        set(gca, "FontSize", 18)
        drawnow limitrate
        hold off
    end
end