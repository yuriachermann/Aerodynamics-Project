function LLT_plot_results(A, flycond, wing, n)

% span 
y = linspace(-wing.b/2, wing.b/2, n);
theta = acos((2*y)/wing.b);
N = length(A);
gamma = zeros(1, n);

for k=1:N
    gamma = gamma + A(k) * sin(k * theta);
end

gamma = gamma * 2 * wing.b * flycond.Uinf;

figure;
plot(y, gamma, 'Linewidth', 2.0, 'Color', [1.0 0.5 0.0]);
xlabel('y [m]'); ylabel('\Gamma [m^3 s^{-1}]'); grid on;

delta = sum((A(2:end)/A(1)).^2);
% quando coloca ponto antes do expoente eleva cada elemento
% se não usasse ponto seria feito multiplicação de matrizes
e = 1 / (1 + delta);

CL = pi * wing.AR * A(1);
CDi = (1/(pi*wing.AR*e)) * CL^2;

% display the results
disp('----------------------------------------');
disp('LLT results');
disp(['  ->Lift coefficient: ' num2str(CL, '%6.4f')])
disp(['  ->Induced drag coefficient: ' num2str(CDi, '%6.4f')])
disp('----------------------------------------');

end

% b = 10; % wing span
% 
% n = 20;
% theta  = linspace (0, pi, n);
% y = -(b/2) * cos(theta);
% gamma1 = sin(theta);
% gamma2 = sin(2*theta);
% gamma3 = sin(3*theta);
% gamma4 = sin(4*theta);
% 
% subplot(4,2,1);plot (y, gamma1);
% xlabel ('y[m]'); ylabel('\Gamma m^3 s^{-1}'); grid on; title('Mode 1');
% 
% subplot(4,2,3);plot (y, gamma2);
% xlabel ('y[m]'); ylabel('\Gamma m^3 s^{-1}'); grid on; title('Mode 2');
% 
% subplot(4,2,5);plot (y, gamma3);
% xlabel ('y[m]'); ylabel('\Gamma m^3 s^{-1}'); grid on; title('Mode 3');
% 
% subplot(4,2,7);plot (y, gamma4);
% xlabel ('y[m]'); ylabel('\Gamma m^3 s^{-1}'); grid on; title('Mode 4');
% 
% subplot(4,2,[2 4]); % plotar várias curvas
% plot (y, gamma1,'k'); hold on; plot (y, gamma2,'r');plot (y, gamma3,'b');plot (y, gamma4,'m');
% legend('Mode 1', 'Mode 2', 'Mode 3', 'Mode 4');
% xlabel ('y[m]'); ylabel('\Gamma m^3 s^{-1}'); grid on; title('All modes');
% 
% A = [1.0 0.25 0.5 0.1]; % for a matrix: B = [1.0 2.0; 3.0 4.0]
% gamma = A(1)*gamma1 + A(2)*gamma2 + A(3)*gamma3 + A(4)*gamma4;
% 
% subplot(4,2,[6 8]);
% plot (y, gamma);
% xlabel ('y[m]'); ylabel('\Gamma m^3 s^{-1}');
% grid on; title('Final circulation');