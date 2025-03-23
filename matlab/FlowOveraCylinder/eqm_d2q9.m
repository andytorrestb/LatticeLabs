function f_eq_d2q9 = eqm_d2q9(Rho,U)
% eqm_d2q9(Rho, U) computes the D2Q9 equilibrium PDF for multiple points.
%   Rho is N×1
%   U   is N×2   (each row is [u, v])
%   f_eq_d2q9 is N×9

% Lattice velocity set: 2×9
vLatt = [ 0,  1,  0, -1,  0,  1, -1, -1,  1; 
          0,  0,  1,  0, -1,  1,  1, -1, -1 ];

% Weights: 1×9
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Number of points
N = size(Rho,1);

% 1) Compute the dot product of each U with each lattice direction.
%    vLatt is (2×9), U' is (2×N), so dotVU will be (9×N).
dotVU = vLatt' * U';     % 9×N

% Transpose to make it N×9 so each row corresponds to one point
dotVU = dotVU';          % N×9

% 2) Square of dotVU, element-wise
dotVU2 = dotVU.^2;       % N×9

% 3) The magnitude-squared of the velocity at each point, Nx1
speed2 = sum(U.^2, 2);   % sum of [u^2 + v^2], Nx1

% 4) Expand Rho, speed2, and w so we can do element-wise operations
Rho   = repmat(Rho,    1, 9);  % N×9
speed2= repmat(speed2, 1, 9);  % N×9
w     = repmat(w,      N, 1);  % N×9

% 5) Form the standard D2Q9 equilibrium:
%        f_eq = w * Rho * [ 1 + 3·(e·u) + 9/2·(e·u)^2 - 3/2·|u|^2 ]
%
%    Where e is each discrete velocity (the columns of vLatt).
f_eq_d2q9 = Rho .* ...
            ( 1 + 3*dotVU + (9/2)*dotVU2 - (3/2)*speed2 ) ...
            .* w;

end
