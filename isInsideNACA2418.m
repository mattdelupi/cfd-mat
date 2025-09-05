function isInside = isInsideNACA2418(x, y, c, le, alpha_deg)

isInside = false;

m = 0.02;
p = 0.40;
t = 0.18 * c;

theta = @(xc) ...
    atan(2*m / p^2 * (p - xc) * (xc <= p && xc >= 0) + ...
    2*m / (1 - p^2) * (p - xc) * (xc > p && xc <= 1));
yt = @(xc) ...
    5 * t * (0.2969 * sqrt(xc) - 0.1260 * xc - 0.3516 * xc.^2 + ...
    0.2843 * xc.^3 - 0.1015 * xc.^4);
yc = @(xc) ...
    m / p^2 * (2 * p * xc - xc.^2) * (xc >= 0 && xc <= p) + ...
        m / (1 - p)^2 * (1 - 2*p + 2*p*xc - xc.^2) * (xc > p && xc <= 1);

xU = @(xc) c * (xc - yt(xc) * sin(theta(xc)));
xL= @(xc) c * (xc + yt(xc) * sin(theta(xc)));

yU = @(xc) c * (yc(xc) + yt(xc) .* cos(theta(xc)));
yL = @(xc) c * (yc(xc) - yt(xc) .* cos(theta(xc)));

Upper = @(xc) Rotate([xU(xc), yU(xc)], alpha_deg) + le(:);
Lower = @(xc) Rotate([xL(xc), yL(xc)], alpha_deg) + le(:);

yUpper = @(xc) [0, 1] * Upper(xc);
yLower = @(xc) [0, 1] * Lower(xc);

xte = le(1) + c * cosd(alpha_deg);

if x >= le(1) && x <= xte
    if y >= yLower((x - le(1)) / (c * cosd(alpha_deg))) && y <= yUpper((x - le(1)) / (c * cosd(alpha_deg)))
        isInside = true;
    end
end

end