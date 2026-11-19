function s7()
% Complete phased-array MATLAB demo combining:
%   1) Original-style spatial field snapshots E(x,y,t)
%   2) Original-style polar array-factor plots
%   3) Chebyshev amplitude tapering
%   4) 3-bit phase quantization
%   5) Linear vs circular array comparison
%   6) Numerical metrics: HPBW, PSLL, pointing error, quantization loss
%   7) Circular-array Chebyshev tapering and phase quantization
%   8) Separate phase-quantization figures for 60 deg and 55 deg steering

clear; clc; close all;
set(0,'defaultFigureColor',[1 1 1]);

%  1. Physical and numerical parameters
freq   = 1e9;          % Hz
c      = 3e8;          % m/s
lambda = c/freq;       % wavelength
T      = 1/freq;       % period
omega  = 2*pi*freq;    % angular frequency
k      = 2*pi/lambda;  % wavenumber

% Angular grid for array-factor plots
theta = linspace(0, 2*pi, 4001);
thetaDeg = theta*180/pi;

% Spatial grid parameters for field snapshots
Ns = 30;               % samples per wavelength
Rmax = 8*lambda;       % field view radius
t0 = 0;                % one time snapshot

%  2. Original-style visualization, like the first demo
%  7-element linear array, d = lambda/2 and d = lambda/4, broadside

N0 = 7;
theta0_original = pi/2;       % broadside for linear array along x-axis
ampUniform7 = ones(N0,1);

[xLin_half,    yLin_half]    = linearArray2D(N0, lambda/2);
[xLin_quarter, yLin_quarter] = linearArray2D(N0, lambda/4);

figure('Name','Original-style field snapshots');
subplot(1,2,1);
plotFieldSnapshotAxes(gca, xLin_half, yLin_half, ampUniform7, theta, ...
    theta0_original, k, omega, t0, lambda, Ns, Rmax, inf, true, ...
    'Antenna separation: d = \lambda/2', ...
    'Relative phase difference: \delta = 0^o');

subplot(1,2,2);
plotFieldSnapshotAxes(gca, xLin_quarter, yLin_quarter, ampUniform7, theta, ...
    theta0_original, k, omega, t0, lambda, Ns, Rmax, inf, true, ...
    'Antenna separation: d = \lambda/4', ...
    'Relative phase difference: \delta = 0^o');

% Original-style polar plot of the array factor
AF_original = arrayFactor2D(xLin_half, yLin_half, ampUniform7, theta, ...
    theta0_original, k, inf);

figure('Name','Original-style polar array factor');
polar(theta, AF_original/max(AF_original));
title('Linear array polar AF, N = 7, d = \lambda/2, \delta = 0^o');

%  3. Extended study parameters
N = 16;                 % number of elements for the extended study
d = lambda/2;           % nominal inter-element spacing
theta0_deg = 60;        % off-broadside steering angle
theta0 = theta0_deg*pi/180;

SLL = 30;               % Chebyshev sidelobe target in dB
B = 3;                  % phase shifter bits
analysisSectorDeg = 180;

if exist('chebwin','file') ~= 2
    error(['chebwin was not found. It requires the Signal Processing Toolbox. ', ...
           'Install it or replace wCheb with another taper.']);
end

%  4. Create geometries
[xLin, yLin]   = linearArray2D(N, d);
[xCirc, yCirc] = circularArray2D(N, d);

%  5. Amplitude weights
wUniform = ones(N,1);

wCheb = chebwin(N, SLL);
wCheb = wCheb / max(wCheb);

%  6. Compute array factors
% Linear array
AF_lin_uniform_ideal = arrayFactor2D(xLin, yLin, wUniform, theta, theta0, k, inf);
AF_lin_cheb_ideal    = arrayFactor2D(xLin, yLin, wCheb,    theta, theta0, k, inf);
AF_lin_uniform_q     = arrayFactor2D(xLin, yLin, wUniform, theta, theta0, k, B);
AF_lin_cheb_q        = arrayFactor2D(xLin, yLin, wCheb,    theta, theta0, k, B);

% Circular array
AF_circ_uniform_ideal = arrayFactor2D(xCirc, yCirc, wUniform, theta, theta0, k, inf);
AF_circ_uniform_q     = arrayFactor2D(xCirc, yCirc, wUniform, theta, theta0, k, B);
AF_circ_cheb_ideal    = arrayFactor2D(xCirc, yCirc, wCheb,    theta, theta0, k, inf);
AF_circ_cheb_q        = arrayFactor2D(xCirc, yCirc, wCheb,    theta, theta0, k, B);

% Extra phase-quantization study for the linear Chebyshev case
B_high = 4;             % higher phase resolution: 16 states, 22.5 deg step
B_low  = 1;             % lower phase resolution: 2 states, 180 deg step
theta0_alt_deg = 55;    % alternative steering angle to make quantization error visible
theta0_alt = theta0_alt_deg*pi/180;

AF_lin_cheb_q_4bit      = arrayFactor2D(xLin, yLin, wCheb, theta, theta0,     k, B_high);
AF_lin_cheb_q_1bit      = arrayFactor2D(xLin, yLin, wCheb, theta, theta0,     k, B_low);
AF_lin_cheb_ideal_55    = arrayFactor2D(xLin, yLin, wCheb, theta, theta0_alt, k, inf);
AF_lin_cheb_q_4bit_55   = arrayFactor2D(xLin, yLin, wCheb, theta, theta0_alt, k, B_high);
AF_lin_cheb_q_3bit_55   = arrayFactor2D(xLin, yLin, wCheb, theta, theta0_alt, k, B);
AF_lin_cheb_q_1bit_55   = arrayFactor2D(xLin, yLin, wCheb, theta, theta0_alt, k, B_low);


%  7. dB plots for analytical comparison
figure('Name','Linear uniform vs Chebyshev');
plot(thetaDeg, normDb(AF_lin_uniform_ideal), 'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_lin_cheb_ideal),    'LineWidth', 1.5);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Linear Array: Uniform vs Chebyshev, \theta_0 = ' num2str(theta0_deg) '^o']);
legend('Uniform', ['Chebyshev, SLL = ' num2str(SLL) ' dB']);

figure('Name','Linear phase quantization');
plot(thetaDeg, normDb(AF_lin_cheb_ideal), 'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_lin_cheb_q), '--', 'LineWidth', 1.5);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Linear Array: Chebyshev Ideal Phase vs ' num2str(B) '-bit Phase, \theta_0 = ' num2str(theta0_deg) '^o']);
legend('Chebyshev ideal phase', ['Chebyshev ' num2str(B) '-bit phase']);

figure('Name','Linear Chebyshev phase bits, theta0 = 60 deg');
plot(thetaDeg, normDb(AF_lin_cheb_ideal),  'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_lin_cheb_q_4bit), '--', 'LineWidth', 1.5);
plot(thetaDeg, normDb(AF_lin_cheb_q),      '-.', 'LineWidth', 1.5);
plot(thetaDeg, normDb(AF_lin_cheb_q_1bit), ':',  'LineWidth', 1.8);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Linear Chebyshev: Ideal vs 4-bit, 3-bit, 1-bit, \theta_0 = ' num2str(theta0_deg) '^o']);
legend( ...
    ['Ideal, \theta_0 = ' num2str(theta0_deg) '^o'], ...
    [num2str(B_high) '-bit, \theta_0 = ' num2str(theta0_deg) '^o'], ...
    [num2str(B)      '-bit, \theta_0 = ' num2str(theta0_deg) '^o'], ...
    [num2str(B_low)  '-bit, \theta_0 = ' num2str(theta0_deg) '^o'] ...
);

figure('Name','Linear Chebyshev phase bits, theta0 = 55 deg');
plot(thetaDeg, normDb(AF_lin_cheb_ideal_55),  'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_lin_cheb_q_4bit_55), '--', 'LineWidth', 1.5);
plot(thetaDeg, normDb(AF_lin_cheb_q_3bit_55), '-.', 'LineWidth', 1.5);
plot(thetaDeg, normDb(AF_lin_cheb_q_1bit_55), ':',  'LineWidth', 1.8);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Linear Chebyshev: Ideal vs 4-bit, 3-bit, 1-bit, \theta_0 = ' num2str(theta0_alt_deg) '^o']);
legend( ...
    ['Ideal, \theta_0 = ' num2str(theta0_alt_deg) '^o'], ...
    [num2str(B_high) '-bit, \theta_0 = ' num2str(theta0_alt_deg) '^o'], ...
    [num2str(B)      '-bit, \theta_0 = ' num2str(theta0_alt_deg) '^o'], ...
    [num2str(B_low)  '-bit, \theta_0 = ' num2str(theta0_alt_deg) '^o'] ...
);

figure('Name','Linear vs circular');
plot(thetaDeg, normDb(AF_lin_uniform_ideal),  'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_circ_uniform_ideal), 'LineWidth', 1.5);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Linear Array vs Circular Array, \theta_0 = ' num2str(theta0_deg) '^o']);
legend('Linear array', 'Circular array');

figure('Name','Circular uniform vs Chebyshev');
plot(thetaDeg, normDb(AF_circ_uniform_ideal), 'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_circ_cheb_ideal),    'LineWidth', 1.5);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Circular Array: Uniform vs Chebyshev, \theta_0 = ' num2str(theta0_deg) '^o']);
legend('Circular uniform', ['Circular Chebyshev, SLL = ' num2str(SLL) ' dB']);

figure('Name','Circular phase quantization');
plot(thetaDeg, normDb(AF_circ_uniform_ideal), 'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_circ_uniform_q), '--', 'LineWidth', 1.5);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Circular Array: Uniform Ideal Phase vs ' num2str(B) '-bit Phase, \theta_0 = ' num2str(theta0_deg) '^o']);
legend('Circular uniform ideal phase', ['Circular uniform ' num2str(B) '-bit phase']);

figure('Name','Circular Chebyshev phase quantization');
plot(thetaDeg, normDb(AF_circ_cheb_ideal), 'LineWidth', 1.5); hold on;
plot(thetaDeg, normDb(AF_circ_cheb_q), '--', 'LineWidth', 1.5);
grid on; ylim([-60 0]); xlim([0 360]);
xlabel('\theta [deg]');
ylabel('Normalized Array Factor [dB]');
title(['Circular Array: Chebyshev Ideal Phase vs ' num2str(B) '-bit Phase, \theta_0 = ' num2str(theta0_deg) '^o']);
legend('Circular Chebyshev ideal phase', ['Circular Chebyshev ' num2str(B) '-bit phase']);

%  8. Extended field snapshots, now including tapering, quantization,
%  and circular geometry
figure('Name','Extended field snapshots');
subplot(1,3,1);
plotFieldSnapshotAxes(gca, xLin, yLin, wCheb, theta, theta0, k, omega, ...
    t0, lambda, Ns, Rmax, B, false, ...
    ['Linear Chebyshev, ' num2str(B) '-bit phase'], ...
    ['\theta_0 = ' num2str(theta0_deg) '^o']);

subplot(1,3,2);
plotFieldSnapshotAxes(gca, xCirc, yCirc, wUniform, theta, theta0, k, omega, ...
    t0, lambda, Ns, Rmax, B, false, ...
    ['Circular uniform, ' num2str(B) '-bit phase'], ...
    ['\theta_0 = ' num2str(theta0_deg) '^o']);

subplot(1,3,3);
plotFieldSnapshotAxes(gca, xCirc, yCirc, wCheb, theta, theta0, k, omega, ...
    t0, lambda, Ns, Rmax, B, false, ...
    ['Circular Chebyshev, ' num2str(B) '-bit phase'], ...
    ['\theta_0 = ' num2str(theta0_deg) '^o']);

%  9. Metrics
m_lin_uniform_ideal  = patternMetrics(theta, AF_lin_uniform_ideal, theta0, analysisSectorDeg);
m_lin_cheb_ideal     = patternMetrics(theta, AF_lin_cheb_ideal,    theta0, analysisSectorDeg);
m_lin_cheb_q         = patternMetrics(theta, AF_lin_cheb_q,        theta0, analysisSectorDeg);
m_lin_cheb_q_4bit    = patternMetrics(theta, AF_lin_cheb_q_4bit,   theta0, analysisSectorDeg);
m_lin_cheb_q_1bit    = patternMetrics(theta, AF_lin_cheb_q_1bit,   theta0, analysisSectorDeg);
m_lin_cheb_ideal_55  = patternMetrics(theta, AF_lin_cheb_ideal_55, theta0_alt, analysisSectorDeg);
m_lin_cheb_q_4bit_55 = patternMetrics(theta, AF_lin_cheb_q_4bit_55, theta0_alt, analysisSectorDeg);
m_lin_cheb_q_3bit_55 = patternMetrics(theta, AF_lin_cheb_q_3bit_55, theta0_alt, analysisSectorDeg);
m_lin_cheb_q_1bit_55 = patternMetrics(theta, AF_lin_cheb_q_1bit_55, theta0_alt, analysisSectorDeg);
m_circ_uniform_ideal = patternMetrics(theta, AF_circ_uniform_ideal, theta0, analysisSectorDeg);
m_circ_uniform_q     = patternMetrics(theta, AF_circ_uniform_q,     theta0, analysisSectorDeg);
m_circ_cheb_ideal    = patternMetrics(theta, AF_circ_cheb_ideal,    theta0, analysisSectorDeg);
m_circ_cheb_q        = patternMetrics(theta, AF_circ_cheb_q,        theta0, analysisSectorDeg);

loss_lin_cheb_q_dB      = quantizationLoss(theta, AF_lin_cheb_ideal, AF_lin_cheb_q, theta0, analysisSectorDeg);
loss_lin_cheb_4bit_dB = quantizationLoss(theta, AF_lin_cheb_ideal, AF_lin_cheb_q_4bit, theta0, analysisSectorDeg);
loss_lin_cheb_1bit_dB = quantizationLoss(theta, AF_lin_cheb_ideal, AF_lin_cheb_q_1bit, theta0, analysisSectorDeg);
loss_lin_cheb_4bit_55_dB = quantizationLoss(theta, AF_lin_cheb_ideal_55, AF_lin_cheb_q_4bit_55, theta0_alt, analysisSectorDeg);
loss_lin_cheb_3bit_55_dB = quantizationLoss(theta, AF_lin_cheb_ideal_55, AF_lin_cheb_q_3bit_55, theta0_alt, analysisSectorDeg);
loss_lin_cheb_1bit_55_dB = quantizationLoss(theta, AF_lin_cheb_ideal_55, AF_lin_cheb_q_1bit_55, theta0_alt, analysisSectorDeg);
loss_circ_q_dB        = quantizationLoss(theta, AF_circ_uniform_ideal, AF_circ_uniform_q, theta0, analysisSectorDeg);
loss_circ_cheb_q_dB   = quantizationLoss(theta, AF_circ_cheb_ideal, AF_circ_cheb_q, theta0, analysisSectorDeg);

fprintf('\n============================================================\n');
fprintf('Final phased-array demo metrics\n');
fprintf('freq = %.3e Hz, lambda = %.3f m\n', freq, lambda);
fprintf('N = %d, d/lambda = %.3f, theta0 = %.1f deg\n', N, d/lambda, theta0_deg);
fprintf('Metrics are evaluated inside +/- %.1f deg around theta0.\n', analysisSectorDeg/2);
fprintf('============================================================\n\n');

printMetrics('Linear uniform ideal', m_lin_uniform_ideal);
printMetrics('Linear Chebyshev ideal', m_lin_cheb_ideal);
printMetrics(['Linear Chebyshev ' num2str(B) '-bit'], m_lin_cheb_q);
fprintf('Linear Chebyshev quantization peak loss: %.3f dB\n\n', loss_lin_cheb_q_dB);

fprintf('Extra linear Chebyshev phase-quantization study, theta0 = %.1f deg\n', theta0_deg);
printMetrics(['Linear Chebyshev ' num2str(B_high) '-bit, theta0 = ' num2str(theta0_deg) ' deg'], m_lin_cheb_q_4bit);
printMetrics(['Linear Chebyshev ' num2str(B)      '-bit, theta0 = ' num2str(theta0_deg) ' deg'], m_lin_cheb_q);
printMetrics(['Linear Chebyshev ' num2str(B_low)  '-bit, theta0 = ' num2str(theta0_deg) ' deg'], m_lin_cheb_q_1bit);
fprintf('Linear Chebyshev 4-bit quantization peak loss at theta0 = %.1f deg: %.3f dB\n', theta0_deg, loss_lin_cheb_4bit_dB);
fprintf('Linear Chebyshev 3-bit quantization peak loss at theta0 = %.1f deg: %.3f dB\n', theta0_deg, loss_lin_cheb_q_dB);
fprintf('Linear Chebyshev 1-bit quantization peak loss at theta0 = %.1f deg: %.3f dB\n\n', theta0_deg, loss_lin_cheb_1bit_dB);

fprintf('Extra linear Chebyshev phase-quantization study, theta0 = %.1f deg\n', theta0_alt_deg);
printMetrics(['Linear Chebyshev ideal, theta0 = ' num2str(theta0_alt_deg) ' deg'], m_lin_cheb_ideal_55);
printMetrics(['Linear Chebyshev ' num2str(B_high) '-bit, theta0 = ' num2str(theta0_alt_deg) ' deg'], m_lin_cheb_q_4bit_55);
printMetrics(['Linear Chebyshev ' num2str(B)      '-bit, theta0 = ' num2str(theta0_alt_deg) ' deg'], m_lin_cheb_q_3bit_55);
printMetrics(['Linear Chebyshev ' num2str(B_low)  '-bit, theta0 = ' num2str(theta0_alt_deg) ' deg'], m_lin_cheb_q_1bit_55);
fprintf('Linear Chebyshev 4-bit quantization peak loss at theta0 = %.1f deg: %.3f dB\n', theta0_alt_deg, loss_lin_cheb_4bit_55_dB);
fprintf('Linear Chebyshev 3-bit quantization peak loss at theta0 = %.1f deg: %.3f dB\n', theta0_alt_deg, loss_lin_cheb_3bit_55_dB);
fprintf('Linear Chebyshev 1-bit quantization peak loss at theta0 = %.1f deg: %.3f dB\n\n', theta0_alt_deg, loss_lin_cheb_1bit_55_dB);

printMetrics('Circular uniform ideal', m_circ_uniform_ideal);
printMetrics(['Circular uniform ' num2str(B) '-bit'], m_circ_uniform_q);
fprintf('Circular uniform quantization peak loss: %.3f dB\n\n', loss_circ_q_dB);

printMetrics('Circular Chebyshev ideal', m_circ_cheb_ideal);
printMetrics(['Circular Chebyshev ' num2str(B) '-bit'], m_circ_cheb_q);
fprintf('Circular Chebyshev quantization peak loss: %.3f dB\n', loss_circ_cheb_q_dB);

fprintf('============================================================\n');

end

% Local functions
function [x, y] = linearArray2D(N, d)
    n = (0:N-1).';
    x = (n - (N-1)/2) * d;
    y = zeros(N,1);
end

function [x, y] = circularArray2D(N, d)
    n = (0:N-1).';

    % Radius chosen so that adjacent chord spacing is approximately d
    R = d / (2*sin(pi/N));

    alpha = 2*pi*n/N;
    x = R*cos(alpha);
    y = R*sin(alpha);
end

function AF = arrayFactor2D(x, y, amp, theta, theta0, k, phaseBits)

    x = x(:);
    y = y(:);
    amp = amp(:);

    phi = steeringPhases(x, y, theta0, k, phaseBits);

    w = amp .* exp(1j*phi);

    AF = zeros(size(theta));

    for it = 1:length(theta)
        spatialPhase = k*(x*cos(theta(it)) + y*sin(theta(it)));
        AF(it) = abs(sum(w .* exp(1j*spatialPhase)));
    end
end

function phi = steeringPhases(x, y, theta0, k, phaseBits)

    phi = -k*(x*cos(theta0) + y*sin(theta0));

    if isfinite(phaseBits)
        phi = quantizePhase(phi, phaseBits);
    end
end

function phiQ = quantizePhase(phi, B)

    L = 2^B;
    deltaPhi = 2*pi/L;

    % Wrap to [-pi, pi]
    phiWrapped = angle(exp(1j*phi));

    % Quantize to nearest available state
    phiQ = deltaPhi * round(phiWrapped/deltaPhi);

    % Wrap again
    phiQ = angle(exp(1j*phiQ));
end

function ydb = normDb(AF)
    AF = abs(AF);
    ydb = 20*log10(AF/max(AF) + eps);
end

function plotFieldSnapshotAxes(ax, xe, ye, amp, theta, theta0, k, omega, ...
    t0, lambda, Ns, Rmax, phaseBits, upperHalfOnly, xLabelText, titleText)

    axes(ax);

    ds = lambda/Ns;
    R = 0:ds:Rmax;

    xg = R.'*cos(theta);
    yg = R.'*sin(theta);

    phi = steeringPhases(xe, ye, theta0, k, phaseBits);

    E = zeros(size(xg));

    for n = 1:length(xe)
        Rn = sqrt((xg - xe(n)).^2 + (yg - ye(n)).^2);
        E = E + amp(n)*cos(omega*t0 - k*Rn + phi(n));
    end

    % Normalize field only for visualization
    E = E / (max(abs(E(:))) + eps);

    normScale = Rmax;

    pcolor(xg/normScale, yg/normScale, E);
    shading interp;
    hold on;

%     % Array factor overlay in normalized polar coordinates
    AF = arrayFactor2D(xe, ye, amp, theta, theta0, k, phaseBits);
    AF = AF / max(AF);

    addPolarGrid(upperHalfOnly);

    plot(AF.*cos(theta), AF.*sin(theta), 'k', 'LineWidth', 2);

    % Element positions
    plot(xe/normScale, ye/normScale, 'k.', 'MarkerSize', 12);

    axis equal;

    if upperHalfOnly
        xlim([-1 1]);
        ylim([0 1]);
    else
        xlim([-1 1]);
        ylim([-1 1]);
    end

    xlabel(xLabelText, 'FontSize', 10);
    title(titleText, 'FontSize', 10, 'Color', [0 0.5 0.5]);

    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    box on;
end

function addPolarGrid(upperHalfOnly)

    th = linspace(0, 2*pi, 720);

    if upperHalfOnly
        thGrid = linspace(0, pi, 361);
    else
        thGrid = th;
    end

    % Radius circles
    radii = [0.5 1.0];
    for r = radii
        plot(r*cos(thGrid), r*sin(thGrid), ':k', 'LineWidth', 0.7);
    end

    % Angle spokes
    if upperHalfOnly
        angleList = 0:30:180;
    else
        angleList = 0:30:330;
    end

    for a = angleList
        ar = a*pi/180;
        plot([0 cos(ar)], [0 sin(ar)], ':k', 'LineWidth', 0.7);
    end

    % Small labels similar to the original plot
    text(0.04, 0.52, '0.5', 'FontSize', 8);
    text(0.04, 1.00, '1', 'FontSize', 8);
end

function loss_dB = quantizationLoss(theta, AFideal, AFq, theta0, sectorDeg)

    rel = wrapToPiLocal(theta - theta0);
    sectorRad = (sectorDeg*pi/180)/2;
    mask = abs(rel) <= sectorRad;

    peakIdeal = max(abs(AFideal(mask)));
    peakQ     = max(abs(AFq(mask)));

    loss_dB = 20*log10(peakQ/peakIdeal + eps);
end

function m = patternMetrics(theta, AF, theta0, sectorDeg)

    AF = abs(AF);
    rel = wrapToPiLocal(theta - theta0);

    sectorRad = (sectorDeg*pi/180)/2;
    mask = abs(rel) <= sectorRad;

    relSect = rel(mask);
    AFSect  = AF(mask);

    [relSort, order] = sort(relSect);
    AFSort = AFSect(order);

    [peakVal, idxPeak] = max(AFSort);
    AFNorm = AFSort / peakVal;
    AFdB = 20*log10(AFNorm + eps);

    peakRel = relSort(idxPeak);
    m.peakAngleDeg = wrapTo360Local((theta0 + peakRel)*180/pi);
    m.pointingErrorDeg = peakRel*180/pi;

    % HPBW: -3 dB crossings
    left = idxPeak;
    while left > 1 && AFdB(left) > -3
        left = left - 1;
    end

    right = idxPeak;
    while right < length(AFdB) && AFdB(right) > -3
        right = right + 1;
    end

    if left == 1 || right == length(AFdB)
        m.HPBWdeg = NaN;
    else
        leftCross = interpCrossing(relSort(left), relSort(left+1), ...
            AFdB(left), AFdB(left+1), -3);
        rightCross = interpCrossing(relSort(right-1), relSort(right), ...
            AFdB(right-1), AFdB(right), -3);
        m.HPBWdeg = (rightCross - leftCross)*180/pi;
    end

    % Approximate first nulls around the main lobe
    leftNull = idxPeak - 1;
    while leftNull > 2
        if AFNorm(leftNull) < AFNorm(leftNull-1) && AFNorm(leftNull) < AFNorm(leftNull+1)
            break;
        end
        leftNull = leftNull - 1;
    end

    rightNull = idxPeak + 1;
    while rightNull < length(AFNorm)-1
        if AFNorm(rightNull) < AFNorm(rightNull-1) && AFNorm(rightNull) < AFNorm(rightNull+1)
            break;
        end
        rightNull = rightNull + 1;
    end

    sllMask = true(size(AFNorm));

    if leftNull > 2 && rightNull < length(AFNorm)-1
        sllMask(leftNull:rightNull) = false;
    else
        if ~isnan(m.HPBWdeg)
            mainHalfWidth = (m.HPBWdeg*pi/180);
            sllMask(abs(relSort - peakRel) <= mainHalfWidth) = false;
        else
            sllMask(max(1,idxPeak-20):min(length(AFNorm),idxPeak+20)) = false;
        end
    end

    if any(sllMask)
        m.PSLLdB = max(AFdB(sllMask));
    else
        m.PSLLdB = NaN;
    end
end

function x = interpCrossing(x1, x2, y1, y2, yTarget)
    if abs(y2-y1) < eps
        x = (x1+x2)/2;
    else
        x = x1 + (yTarget-y1)*(x2-x1)/(y2-y1);
    end
end

function printMetrics(name, m)
    fprintf('%s\n', name);
    fprintf('  Peak angle      : %.2f deg\n', m.peakAngleDeg);
    fprintf('  Pointing error  : %.2f deg\n', m.pointingErrorDeg);
    fprintf('  HPBW            : %.2f deg\n', m.HPBWdeg);
    fprintf('  PSLL            : %.2f dB\n\n', m.PSLLdB);
end

function y = wrapToPiLocal(x)
    y = mod(x + pi, 2*pi) - pi;
end

function y = wrapTo360Local(xDeg)
    y = mod(xDeg, 360);
end
