function dx = cells(t, x, n_par)
% cells - Model cell expansion kinetics in a bioreactor
%
% Inputs:
%   t     - Time (not used here, needed for ODE solvers)
%   x     - State vector:
%           x(1) = Viable cells (cells/mL)
%           x(2) = Total cells (cells/mL)
%           x(3) = Glucose concentration (g/L)
%   n_par - Parameters vector:
%           n_par(1) = MaxGrowth rate (1/h)
%           n_par(2) = Monod constant for glucose (kglc, g/L)
%           n_par(3) = Maintenance coefficient (mglc, g/(cells*h))
%           n_par(4) = Yield coefficient Yx/glc (cells/g)
%
% Outputs:
%   dx - Rate of change of state vector [dx1; dx2; dx3]

% ------------------------------
% Parameters
% ------------------------------
MaxGrowth = n_par(1);   % Maximum specific growth rate
kglc      = n_par(2);   % Glucose saturation constant
mglc      = n_par(3);   % Maintenance coefficient
Yxglc     = n_par(4);   % Yield coefficient (cells per g glucose)

% ------------------------------
% Variables
% ------------------------------
ViableCell = x(1);
TotalCell  = x(2);
Glc        = x(3);

% ------------------------------
% Bioreactor settings
% ------------------------------
V       = 3;    % Reactor volume (L)
Flowin  = 0;    % Inlet flow (L/h)
Flowout = 0;    % Outlet flow (L/h)
GLCin   = 0;    % Inlet glucose concentration (g/L)

% ------------------------------
% Growth and substrate consumption
% ------------------------------
% Monod kinetics for glucose-limited growth
flim = Glc / (kglc + Glc);  

% Specific growth rate
GrowthRate = MaxGrowth * flim;

% Glucose consumption rate per cell (includes maintenance)
Qglc = GrowthRate / Yxglc + mglc;

% ------------------------------
% Differential equations
% ------------------------------
dx1 = GrowthRate * ViableCell - (Flowout/V) * ViableCell;   % Viable cells
dx2 = GrowthRate * ViableCell - (Flowout/V) * TotalCell;    % Total cells
dx3 = -Qglc * ViableCell + (Flowin/V) * GLCin - (Flowout/V) * Glc;  % Glucose

% Output as column vector
dx = [dx1; dx2; dx3];

end
