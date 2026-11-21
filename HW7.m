%% Cauchy (small) strain of a triangular element (nodes A, B, C)
% With rigid body rotation detection/removal and a full numerical summary.
%
%  1. Ask for undeformed coordinates   2. Ask for deformed coordinates
%  3. Cauchy strain definition         4-5. Displacements & vectors
%  6-9. Displacement field, solve for coefficients
%  Detect rigid body rotation (F = R*U) and undo it if > threshold
%  10. Final strain matrix             11. Rotate coordinate axes
%  12. F = R*U report                  + RESULTS SUMMARY at the end
clear; clc; close all;

rot_tol_deg = 1;   % rotations larger than this (degrees) are removed

%% 1. Undeformed coordinates (re-asks if the points are collinear)
fprintf('--- UNDEFORMED coordinates ---\n');
while true
    A = input('Enter A = [X_A Y_A]: ');
    B = input('Enter B = [X_B Y_B]: ');
    C = input('Enter C = [X_C Y_C]: ');
    area0 = 0.5*((B(1)-A(1))*(C(2)-A(2)) - (C(1)-A(1))*(B(2)-A(2)));
    if abs(area0) > 1e-12, break; end
    fprintf('\n*** A, B, C lie on a straight line (zero area). Enter 3 non-collinear points. ***\n\n');
end

%% 2. Deformed coordinates (matrix form)
fprintf('\n--- DEFORMED coordinates ---\n');
a = input('Enter a = [x_A y_A]: ');
b = input('Enter b = [x_B y_B]: ');
c = input('Enter c = [x_C y_C]: ');

P0 = [A; B; C];   % undeformed, one node per row
P1 = [a; b; c];   % deformed,   one node per row

%% 3. Definition (Cauchy strain)
%  eps = 1/2 * ( grad(u) + grad(u)' )
%  exx = du_x/dx,  eyy = du_y/dy,  exy = 1/2*(du_x/dy + du_y/dx)

%% 4-5. Displacements and displacement vectors (as measured)
dA = (a - A)';   dB = (b - B)';   dC = (c - C)';

%% 6-9. Displacement field  u = a1 + a2*X + a3*Y ,  v = b1 + b2*X + b3*Y
[acoef_raw, bcoef_raw, H_raw] = fit_field(P0, P1);
E_raw = 0.5*(H_raw + H_raw');    % strain WITHOUT removing rotation

%% Detect rigid body rotation  (polar decomposition  F = R*U)
F = eye(2) + H_raw;              % deformation gradient of the triangle
if det(F) <= 0
    warning('det(F) <= 0: triangle is inverted/degenerate, rotation not reliable.');
end
U = sqrtm(F' * F);               % right stretch tensor
R = F / U;                       % rotation tensor
theta_rb = atan2d(R(2,1), R(1,1));   % rigid body rotation angle (deg, CCW +)

if abs(theta_rb) > rot_tol_deg
    xc    = mean(P1, 1);                    % rotate back about deformed centroid
    P1_nr = ((R' * (P1 - xc)')') + xc;      % de-rotated deformed coordinates
    derotated = true;
else
    P1_nr = P1;
    derotated = false;
end

%% Re-fit the displacement field with the rotation removed
[acoef, bcoef, H] = fit_field(P0, P1_nr);

%% 7. Partial derivatives -> strain components
dudx = acoef(2);   % a2
dudy = acoef(3);   % a3
dvdx = bcoef(2);   % b2
dvdy = bcoef(3);   % b3

exx = dudx;
eyy = dvdy;
exy = 0.5*(dudy + dvdx);

%% 10. Final strain
E = [exx exy;
     exy eyy];

% Principal strains and directions
[V, D] = eig(E);
[e_princ, idx] = sort(diag(D), 'descend');
V = V(:, idx);
theta_p = atan2d(V(2,1), V(1,1));            % direction of e1 (deg from +X)
gamma_max = e_princ(1) - e_princ(2);         % max in-plane engineering shear

%% 11. Rotated coordinate axes:  E' = Q E Q'
theta_deg = input('\nEnter rotation angle of axes theta (degrees): ');
th = deg2rad(theta_deg);
Q  = [ cos(th)  sin(th);
      -sin(th)  cos(th)];
Eprime = Q * E * Q';

%% Normal components (Cauchy formula in matrix notation: t = E * n)
n = input('Enter direction vector n = [nx ny] for normal components: ');
n = n(:) / norm(n);
t = E * n;
e_normal = n' * E * n;
e_shear  = norm(t - e_normal*n);

%% Save everything to one struct in the workspace
results.undeformed        = P0;
results.deformed          = P1;
results.deformed_derot    = P1_nr;
results.displacements     = [dA dB dC]';
results.a_coef            = acoef;
results.b_coef            = bcoef;
results.F                 = F;
results.R                 = R;
results.U                 = U;
results.rotation_deg      = theta_rb;
results.rotation_removed  = derotated;
results.E_raw             = E_raw;
results.E                 = E;
results.E_rotated_axes    = Eprime;
results.principal_strains = e_princ;
results.principal_dir_deg = theta_p;
results.normal_strain_n   = e_normal;
results.shear_strain_n    = e_shear;

%% ===================== NUMERICAL RESULTS =====================
fmt = '%12.6f';
clc;
fprintf('==============================================================\n');
fprintf('                    NUMERICAL RESULTS\n');
fprintf('==============================================================\n');

fprintf('\n--- Node coordinates ---\n');
fprintf('%-6s %12s %12s %12s %12s\n', 'Node', 'X', 'Y', 'x', 'y');
names = 'ABC';
for k = 1:3
    fprintf('%-6s %12.6f %12.6f %12.6f %12.6f\n', names(k), P0(k,1), P0(k,2), P1(k,1), P1(k,2));
end
fprintf('Triangle area (undeformed) = %.6f\n', abs(area0));

fprintf('\n--- Displacements (measured) ---\n');
fprintf('%-6s %12s %12s\n', 'Node', 'u', 'v');
D3 = [dA dB dC]';
for k = 1:3
    fprintf('%-6s %12.6f %12.6f\n', names(k), D3(k,1), D3(k,2));
end

fprintf('\n--- Rigid body rotation ---\n');
fprintf('Detected rotation = %.6f deg (tolerance %.2f deg)\n', theta_rb, rot_tol_deg);
if derotated
    fprintf('Rotation EXCEEDED tolerance -> it was removed before computing strain.\n');
else
    fprintf('Rotation within tolerance -> no correction applied.\n');
end
fprintf('\nDe-rotated deformed coordinates:\n');
fprintf('%-6s %12s %12s\n', 'Node', 'x''', 'y''');
for k = 1:3
    fprintf('%-6s %12.6f %12.6f\n', names(k), P1_nr(k,1), P1_nr(k,2));
end

fprintf('\n--- Displacement field coefficients (rotation removed) ---\n');
fprintf('u = a1 + a2 X + a3 Y :  a1 = %9.6f   a2 = %9.6f   a3 = %9.6f\n', acoef);
fprintf('v = b1 + b2 X + b3 Y :  b1 = %9.6f   b2 = %9.6f   b3 = %9.6f\n', bcoef);

fprintf('\n--- STRAIN COMPONENTS (rotation removed) ---\n');
fprintf('exx = %s\n', num2str(exx, '%.6g'));
fprintf('eyy = %s\n', num2str(eyy, '%.6g'));
fprintf('exy = %s   (engineering shear gamma_xy = %.6g)\n', num2str(exy, '%.6g'), 2*exy);
fprintf('\nStrain tensor E =\n');
print_matrix(E, fmt);

if derotated
    fprintf('\nFor comparison, strain tensor WITHOUT removing rotation (contains error):\n');
    print_matrix(E_raw, fmt);
    fprintf('Difference (raw - corrected):\n');
    print_matrix(E_raw - E, fmt);
end

fprintf('\n--- Principal strains ---\n');
fprintf('e1 = %.6g   (direction %.4f deg from +X)\n', e_princ(1), theta_p);
fprintf('e2 = %.6g\n', e_princ(2));
fprintf('Max in-plane shear (engineering) = %.6g\n', gamma_max);

fprintf('\n--- Strain in rotated axes (theta = %g deg) ---\n', theta_deg);
fprintf('Q =\n');          print_matrix(Q, fmt);
fprintf('E'' = Q E Q^T =\n'); print_matrix(Eprime, fmt);

fprintf('\n--- Normal components for n = [%.4f %.4f] ---\n', n);
fprintf('t = E*n          = [%.6g ; %.6g]\n', t);
fprintf('Normal strain    = %.6g\n', e_normal);
fprintf('Shear component  = %.6g\n', e_shear);

fprintf('\n--- Polar decomposition F = R*U ---\n');
fprintf('F =\n'); print_matrix(F, fmt);
fprintf('R (rigid rotation) =\n'); print_matrix(R, fmt);
fprintf('U (pure stretch) =\n');   print_matrix(U, fmt);
fprintf('Check |F - R*U| = %.3g\n', norm(F - R*U));

fprintf('\n(All results are also stored in the workspace variable "results".)\n');
fprintf('==============================================================\n');

%% Plot (with strain values written on the figure)
figure('Position', [100 100 800 600]); hold on; axis equal; grid on;
h1 = fill(P0(:,1),    P0(:,2),    [0.8 0.8 1], 'FaceAlpha', 0.5);
h2 = fill(P1(:,1),    P1(:,2),    [1 0.8 0.8], 'FaceAlpha', 0.5);
h3 = fill(P1_nr(:,1), P1_nr(:,2), [0.8 1 0.8], 'FaceAlpha', 0.5);
for k = 1:3
    text(P0(k,1), P0(k,2), ['  ' names(k)]);
    text(P1(k,1), P1(k,2), ['  ' lower(names(k))]);
end
legend([h1 h2 h3], 'Undeformed', 'Deformed (measured)', 'Deformed (rotation removed)', ...
       'Location', 'best');
title(sprintf('\\epsilon_{xx}=%.4g   \\epsilon_{yy}=%.4g   \\epsilon_{xy}=%.4g   (rotation = %.3g^\\circ)', ...
      exx, eyy, exy, theta_rb));
xlabel('X'); ylabel('Y');

%% ---------- Local functions (must stay at end of script) ----------
function [acoef, bcoef, H] = fit_field(P0, P1)
% Fits u = a1 + a2*X + a3*Y and v = b1 + b2*X + b3*Y through 3 nodes.
% P0, P1: 3x2 matrices (one node per row), undeformed and deformed.
    M = [ones(3,1) P0(:,1) P0(:,2)];
    if abs(det(M)) < 1e-12
        error('Undeformed points are collinear (zero area triangle).');
    end
    u = P1(:,1) - P0(:,1);
    v = P1(:,2) - P0(:,2);
    acoef = M \ u;      % [a1; a2; a3]
    bcoef = M \ v;      % [b1; b2; b3]
    H = [acoef(2) acoef(3);
         bcoef(2) bcoef(3)];   % displacement gradient
end

function print_matrix(Mx, fmt)
% Prints a matrix with fixed-width columns.
    for i = 1:size(Mx,1)
        fprintf('   [');
        fprintf(fmt, Mx(i,:));
        fprintf(' ]\n');
    end
end
