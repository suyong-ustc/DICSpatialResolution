function [u, e] = DICStepDisplacementStrain(x, x0, m, N, weight)
% DICStepDisplacementStrain
%   此处显示详细说明

% 无量纲量
z = (x - x0) / m;

% 估计位移场
if weight == 0  % 均匀权重

    if N == 0 || N == 1
        u = 0.5 * (1 + z);
    elseif N == 2 || N == 3
        u = 0.5 + 9/8 * z - 5/8 * z.^3;
    end
    u(z < -1) = 0;
    u(z > 1) = 1;
elseif weight == 1  % 高斯权重

    if  N == 0 || N == 1
        u = 0.5 * (1 + erf(z));
    elseif N == 2 || N == 3
        u = 0.5 + 0.5 * erf(z) + 0.5/sqrt(pi) * z .* exp(-z.*z);
    end

end

% 估计应变场
if weight == 0  % 均匀权重

    if N == 0 || N == 1
        e = 1/(2*m) * ones(size(z));
    elseif N == 2 || N == 3
        e = 3/(8*m) * (3 - 5*z.^2);
    end
    e(abs(z) > 1) = 0;

elseif weight == 1  % 高斯权重

    if  N == 0 || N == 1
        e = 1/sqrt(pi) / m * exp(-z.*z);
    elseif N == 2 || N == 3
        e = 1/sqrt(pi) / m * exp(-z.*z) .* (1.5 - z.*z);
    end

end
   
end