function results = s9()
% Complete phased-array MATLAB demo combining:
%   1) Original-style spatial field snapshots E(x,y,t)
%   2) Original-style polar array-factor plots
%   3) Chebyshev amplitude tapering
%   4) 3-bit reference comparison and 1-5 bit phase-quantization study
%   5) Linear vs circular array comparison
%   6) Numerical metrics: HPBW, PSLL, pointing error, quantization loss
%   7) Circular-array Chebyshev tapering and phase quantization
%   8) Side-by-side 55/60 deg phase-quantization plots, metrics, and table
% All studies use the SAME array-centered phase reference as s7.m.

clear; clc; close all;
set(0,'defaultFigureColor',[1 1 1]);

%  1. Physical and numerical parameters
freq   = 28e9;          % Hz
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
B = 3;                  % 3-bit baseline for the original comparisons
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

% Extended phase-quantization study: linear Chebyshev, 1-5 bits.
% IMPORTANT: Keep the array-centered phase reference used in s7.m.
% No first-element phase subtraction is applied at any resolution.
% Evaluate over the FRONT half-plane to avoid the mirrored main beam of
% a linear array (the full 0-360 deg pattern contains front/back ambiguity).
% Use a fine 0.01-deg angular sampling grid for accurate numerical metrics.
bitValues = 1:5;
steeringAnglesDeg = [55 60];
thetaStudy = linspace(0, pi, 18001);
thetaStudyDeg = thetaStudy*180/pi;
nAngles = numel(steeringAnglesDeg);
nBits = numel(bitValues);

AFstudyIdeal = cell(nAngles,1);
AFstudyQuant = cell(nAngles,nBits);
metricsIdeal = cell(nAngles,1);
metricsQuant = cell(nAngles,nBits);
lossPeak_dB = zeros(nAngles,nBits);
lossTarget_dB = zeros(nAngles,nBits);
psll_dB = zeros(nAngles,nBits);
pointingError_deg = zeros(nAngles,nBits);
hpbw_deg = zeros(nAngles,nBits);

for ia = 1:nAngles
    steer = steeringAnglesDeg(ia)*pi/180;
    AFstudyIdeal{ia} = arrayFactor2D(xLin, yLin, wCheb, thetaStudy, steer, k, inf);
    % sectorDeg = 360 includes the entire sampled front half-plane.
    metricsIdeal{ia} = patternMetrics(thetaStudy, AFstudyIdeal{ia}, steer, 360);

    [~, idxTarget] = min(abs(thetaStudy - steer));
    idealPeak = max(AFstudyIdeal{ia});

    for ib = 1:nBits
        bits = bitValues(ib);
        AFstudyQuant{ia,ib} = arrayFactor2D(xLin, yLin, wCheb, ...
                                             thetaStudy, steer, k, bits);
        metricsQuant{ia,ib} = patternMetrics(thetaStudy, ...
                                             AFstudyQuant{ia,ib}, steer, 360);

        % dB loss of the actual pattern peak, relative to the ideal peak.
        lossPeak_dB(ia,ib) = 20*log10(idealPeak / ...
                                      max(AFstudyQuant{ia,ib}));
        % dB loss exactly at the REQUESTED steering angle.
        lossTarget_dB(ia,ib) = 20*log10(idealPeak / ...
                                        AFstudyQuant{ia,ib}(idxTarget));
        psll_dB(ia,ib) = metricsQuant{ia,ib}.PSLLdB;
        pointingError_deg(ia,ib) = metricsQuant{ia,ib}.pointingErrorDeg;
        hpbw_deg(ia,ib) = metricsQuant{ia,ib}.HPBWdeg;
    end
end

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

% Both angles in one figure: fixed ideal-peak reference is preferred
% because it shows gain loss AND distortion simultaneously. No duplicate
% self-normalized figure (which hides gain loss) is generated.
lineTypes = {'--','-.',':','--','-'};
lineColors = lines(nBits);
legendText = [{'Ideal'}, arrayfun(@(b) sprintf('%d-bit',b), ...
                                   bitValues, 'UniformOutput', false)];
figure('Name','Linear Chebyshev: 55 and 60 deg (same centered reference)');
for ia = 1:nAngles
    subplot(1,2,ia);
    idealPeak = max(AFstudyIdeal{ia});
    plot(thetaStudyDeg, 20*log10(AFstudyIdeal{ia}/idealPeak + eps), ...
         'k-', 'LineWidth', 2);
    hold on;
    for ib = 1:nBits
        plot(thetaStudyDeg, ...
             20*log10(AFstudyQuant{ia,ib}/idealPeak + eps), ...
             'LineStyle', lineTypes{ib}, 'Color', lineColors(ib,:), ...
             'LineWidth', 1.5);
    end
    grid on; xlim([0 180]); ylim([-60 2]);
    xlabel('Angle [deg]');
    ylabel('Array factor relative to ideal peak [dB]');
    title(sprintf('Linear Chebyshev, steering = %d deg', ...
                  steeringAnglesDeg(ia)));
    legend(legendText, 'Location', 'southwest');
end

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

%  9. Metrics for the original comparisons (3-bit baseline)
m_lin_uniform_ideal  = patternMetrics(theta, AF_lin_uniform_ideal, theta0, analysisSectorDeg);
m_lin_cheb_ideal     = patternMetrics(theta, AF_lin_cheb_ideal,    theta0, analysisSectorDeg);
m_lin_cheb_q         = patternMetrics(theta, AF_lin_cheb_q,        theta0, analysisSectorDeg);
m_circ_uniform_ideal = patternMetrics(theta, AF_circ_uniform_ideal, theta0, analysisSectorDeg);
m_circ_uniform_q     = patternMetrics(theta, AF_circ_uniform_q,     theta0, analysisSectorDeg);
m_circ_cheb_ideal    = patternMetrics(theta, AF_circ_cheb_ideal,    theta0, analysisSectorDeg);
m_circ_cheb_q        = patternMetrics(theta, AF_circ_cheb_q,        theta0, analysisSectorDeg);

loss_lin_cheb_q_dB  = quantizationLoss(theta, AF_lin_cheb_ideal, ...
                                      AF_lin_cheb_q, theta0, analysisSectorDeg);
loss_circ_q_dB      = quantizationLoss(theta, AF_circ_uniform_ideal, ...
                                      AF_circ_uniform_q, theta0, analysisSectorDeg);
loss_circ_cheb_q_dB = quantizationLoss(theta, AF_circ_cheb_ideal, ...
                                      AF_circ_cheb_q, theta0, analysisSectorDeg);

fprintf('\n============================================================\n');
fprintf('Final phased-array demo metrics\n');
fprintf('freq = %.3e Hz, lambda = %.3f m\n', freq, lambda);
fprintf('N = %d, d/lambda = %.3f, baseline steering = %.1f deg\n', ...
        N, d/lambda, theta0_deg);
fprintf('Original 3-bit metrics use a +/- %.1f deg angular sector.\n', ...
        analysisSectorDeg/2);
fprintf('============================================================\n\n');

printMetrics('Linear uniform ideal', m_lin_uniform_ideal);
printMetrics('Linear Chebyshev ideal', m_lin_cheb_ideal);
printMetrics(sprintf('Linear Chebyshev %d-bit', B), m_lin_cheb_q);
fprintf('Linear Chebyshev %d-bit peak gain loss: %.3f dB\n\n', ...
        B, loss_lin_cheb_q_dB);

% 10. Comparison table, including ideal reference for both angles.
% 'PeakGainLoss' compares the peaks of the actual and ideal patterns;
% 'TargetGainLoss' compares their response at the desired steering angle.
nRows = nAngles * (nBits + 1);
steeringCol = zeros(nRows,1);
caseCol = cell(nRows,1);
peakLossCol = zeros(nRows,1);
targetLossCol = zeros(nRows,1);
psllCol = zeros(nRows,1);
errorCol = zeros(nRows,1);
hpbwCol = zeros(nRows,1);

ir = 0;
for ia = 1:nAngles
    ir = ir + 1;
    steeringCol(ir) = steeringAnglesDeg(ia);
    caseCol{ir} = 'Ideal';
    peakLossCol(ir) = 0;
    targetLossCol(ir) = 0;
    psllCol(ir) = metricsIdeal{ia}.PSLLdB;
    errorCol(ir) = metricsIdeal{ia}.pointingErrorDeg;
    hpbwCol(ir) = metricsIdeal{ia}.HPBWdeg;

    for ib = 1:nBits
        ir = ir + 1;
        steeringCol(ir) = steeringAnglesDeg(ia);
        caseCol{ir} = sprintf('%d-bit', bitValues(ib));
        peakLossCol(ir) = lossPeak_dB(ia,ib);
        targetLossCol(ir) = lossTarget_dB(ia,ib);
        psllCol(ir) = psll_dB(ia,ib);
        errorCol(ir) = pointingError_deg(ia,ib);
        hpbwCol(ir) = hpbw_deg(ia,ib);
    end
end

results = table(steeringCol, caseCol, peakLossCol, targetLossCol, ...
                psllCol, errorCol, hpbwCol, ...
                'VariableNames', {'Steering_deg', 'Resolution', ...
                'PeakGainLoss_dB', 'TargetGainLoss_dB', 'PSLL_dB', ...
                'PointingError_deg', 'HPBW_deg'});

fprintf('============================================================\n');
fprintf('LINEAR CHEBYSHEV: 1-5 BITS AT 55 AND 60 DEG\n');
fprintf('Metrics use the front half-plane 0-180 deg, 0.01 deg grid.\n');
fprintf('Gain loss is POSITIVE for degradation, PSLL is negative.\n');
fprintf('Signed pointing error = nearest equal-height main peak - target angle.\n');
fprintf('NOTE: 1-bit phase can produce a symmetric equal-height beam (PSLL near 0 dB).\n');
fprintf('Phase reference and quantizer match s7.m (array centered, MATLAB round).\n');
fprintf('============================================================\n');
disp(results);

figure('Name','Quantization metrics at 55 and 60 deg (1-5 bits)');
angleColors = lines(nAngles);
subplot(1,3,1);
hold on;
for ia = 1:nAngles
    plot(bitValues, lossPeak_dB(ia,:), '-o', 'LineWidth', 1.8, ...
         'Color', angleColors(ia,:), 'MarkerSize', 6);
end
grid on; xlim([0.8 5.2]); xticks(bitValues);
xlabel('Phase resolution [bits]'); ylabel('Peak gain loss [dB]');
title('Gain loss');
legend({'55 deg','60 deg'}, 'Location','northeast');

subplot(1,3,2);
hold on;
for ia = 1:nAngles
    plot(bitValues, psll_dB(ia,:), '-o', 'LineWidth', 1.8, ...
         'Color', angleColors(ia,:), 'MarkerSize', 6);
end
grid on; xlim([0.8 5.2]); xticks(bitValues);
xlabel('Phase resolution [bits]'); ylabel('PSLL [dB]');
title('Peak sidelobe level');
legend({'55 deg','60 deg'}, 'Location','northeast');

subplot(1,3,3);
hold on;
for ia = 1:nAngles
    plot(bitValues, abs(pointingError_deg(ia,:)), '-o', ...
         'LineWidth', 1.8, 'Color', angleColors(ia,:), ...
         'MarkerSize', 6);
end
grid on; xlim([0.8 5.2]); xticks(bitValues);
xlabel('Phase resolution [bits]'); ylabel('|Pointing error| [deg]');
title('Absolute pointing error');
legend({'55 deg','60 deg'}, 'Location','northeast');

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

    % Same array-centered reference as s7.m.
    % No subtraction of phi(1): all phases are quantized as in s7.m.
    phi = -k*(x*cos(theta0) + y*sin(theta0));

    if isfinite(phaseBits)
        phi = quantizePhase(phi, phaseBits);
    end
end

function phiQ = quantizePhase(phi, B)

    L = 2^B;
    deltaPhi = 2*pi/L;

    % Same quantizer as s7.m, including MATLAB round() tie behavior.
    phiWrapped = angle(exp(1j*phi));
    phiQ = deltaPhi * round(phiWrapped/deltaPhi);
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

    % Positive dB value = gain reduction due to phase quantization.
    loss_dB = 20*log10(peakIdeal / (peakQ + eps));
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
    % For tied (or numerically identical) maxima, select the beam closest
    % to the requested direction. This is essential for 1-bit real weights,
    % which create symmetric equal-amplitude beams in this geometry.
    tiedPeaks = find(AFSort >= peakVal*(1-1e-10));
    if numel(tiedPeaks) > 1
        [~, closestIdx] = min(abs(relSort(tiedPeaks)));
        idxPeak = tiedPeaks(closestIdx);
    end
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
