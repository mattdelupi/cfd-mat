function dxdt = PathlinesRHS(~, xP, x, y, U, V, Dx, Dy, Nx, Ny)

i = floor(xP(1) / Dx) + 1;
j = floor(xP(2) / Dy) + 1;

xNode = [x(min(Nx, i)); y(j)];

csix = (xP(1) - xNode(1)) / Dx;
csiy = (xP(2) - xNode(2)) / Dy;

if i+1 <= Nx && j+1 <= Ny && i > 0 && j > 0
    u = [1-csix, csix] * U([i, i+1], [j, j+1]) * [1-csiy; csiy];
    v = [1-csix, csix] * V([i, i+1], [j, j+1]) * [1-csiy; csiy];
else
    u = 0.05;
    v = 0;
end

dxdt = [u; v];

end