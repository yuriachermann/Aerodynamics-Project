function A = LLT_solve(wing, flycond, N)
% LLT_solve computes the wing circulation distribution
% The function returns the Fourier series amplitudes
% 
% ARGUMENTS:
% 
%       wing: wing struct contains...
%       flycond: ...
% 
% 
% EXAMPLE: A = LLT_solve(wing, flycond, N)
% 
% See also: LLT_plot_results

theta = linspace(0+pi/(20*N), pi-pi/(20*N), N);
%discretization from 0+qsi to pi-qsi -> here qsi=pi/(20*N)

for i=1:N
    for j = 1:N
        y = -wing.b * cos(theta(j)) / 2;
        c = interp1(wing.stations, wing.chords, y);
        nu = c * wing.airfoil.clalpha / (4 * wing.b);
        M(i,j) = sin(i * theta(j)) .* (nu * i + sin(theta(j)));
        b(j,1) = nu * (wing.alpha - wing.airfoil.alpha_L_0) * sin(theta(j)) * pi / 180;
    end
end

A = inv(M') * b;

end