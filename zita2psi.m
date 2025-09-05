function [PSI, zita] = zita2psi(zita, ZITA, PSI, Lap, G, Nx, Ny, hq, ...
                        PsiS, PsiN, PsiW, PsiE, PsiCyl)

zita(G(G > 0)) = ZITA(G > 0);

zita(G(2:end-1, 2)) = zita(G(2:end-1, 2)) - PsiS(2:end-1) / hq;
zita(G(2:end-1, end-1)) = zita(G(2:end-1, end-1)) - PsiN(2:end-1) / hq;
zita(G(2, 2:end-1)) = zita(G(2, 2:end-1)) - PsiW(2:end-1)' / hq;
zita(G(end-1, 2:end-1)) = zita(G(end-1, 2:end-1)) - PsiE(2:end-1)' / hq;

i = 2 : Nx-1;
j = 2 : Ny-1;
D = G(i, j);
zita(D(and(G(i, j) > 0, G(i, j+1) == -1))) = ...
    zita(D(and(G(i, j) > 0, G(i, j+1) == -1))) ...
    - PsiCyl / hq;
zita(D(and(G(i, j) > 0, G(i, j-1) == -1))) = ...
    zita(D(and(G(i, j) > 0, G(i, j-1) == -1))) ...
    - PsiCyl / hq;
zita(D(and(G(i, j) > 0, G(i+1, j) == -1))) = ...
    zita(D(and(G(i, j) > 0, G(i+1, j) == -1))) ...
    - PsiCyl / hq;
zita(D(and(G(i, j) > 0, G(i-1, j) == -1))) = ...
    zita(D(and(G(i, j) > 0, G(i-1, j) == -1))) ...
    - PsiCyl / hq;

psi = Lap \ zita;
PSI(G > 0) = psi;

end