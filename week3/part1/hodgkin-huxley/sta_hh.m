% Specify the properties of ion channels embedded in the membrane
Na = make_channel(50, 1.2e3);
K = make_channel(-77, 0.36e3);
leak = make_channel(-54.387, 3);

% Specify the membrane parameters
Cm = 10;  % nF/mm^2;

% Simulation time
dt = 0.05;                  % integration time step (s)
T = 500000;% [[ YOUR CODE HERE ]]  % simulation time (s)

% noise
I_sigma = 80% [[ YOUR CODE HERE ]]

% generate the current pulse
time = 0:dt:T;
current = I_sigma * randn(size(time));

% initialize HH variables
voltage = zeros(length(time), 1); voltage(1) = -64.996;
m = zeros(length(time), 1); m(1) = 0.052955;
n = zeros(length(time), 1); n(1) = 0.31773;
h = zeros(length(time), 1); h(1) = 0.59599;

% run the simulation
for t = 2:length(time)

  % update on progress
  if mod(t, 1000/dt) == 0
    fprintf('%f%% done\n', 100 * t / length(time));
  end

  % update gating variables
  [alpha, beta] = gates(voltage(t-1));
  m(t) = m(t-1) + dt * (alpha.m * (1-m(t-1)) - beta.m * m(t-1));
  n(t) = n(t-1) + dt * (alpha.n * (1-n(t-1)) - beta.n * n(t-1));
  h(t) = h(t-1) + dt * (alpha.h * (1-h(t-1)) - beta.h * h(t-1));

  % update voltage
  total_current = current(t-1) + ...
                  I_sigma * randn + ...
                  leak.gmax * (leak.E - voltage(t-1)) + ...
                  Na.gmax * m(t-1)^3 * h(t-1) * (Na.E - voltage(t-1)) + ...
                  K.gmax * n(t-1)^4 * (K.E - voltage(t-1));

  voltage(t) = voltage(t-1) + (dt / Cm) * total_current;

end
%%

figure(1); clf;
subplot(3, 1, 1)
plot(time, current, 'k-');
xlim([time(1) time(end)]);
ylabel('Current (nA)');
set(gca, 'xtick', []);
subplot(3, 1, [2,3])
plot(time, voltage);
xlim([time(1) time(end)]);
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title("Simulation of the HH neuron")
%%
% compute the STA
tau = 200;
spike_times = peakdet(voltage, 30);% [[ YOUR CODE HERE ]]     % identify spikes using peakdet

hold on
plot(time(spike_times(:,1)), spike_times(:,2), '*')
hold off

%%
nspk = length(spike_times);
ste = zeros(3 * tau, nspk);
for idx = 1:nspk
  mi = spike_times(idx, 1);
  if mi >= 2*tau & mi+tau <= length(current)
    ste(:, idx) = current(mi-2*tau+1:mi+tau);
  end
end
sta = mean(ste, 2);
lags = (-2*tau+1:tau) * dt;

figure(2);
subplot(2,2,[1,3])
plot(lags(lags < -0.1), sta(lags < -0.1));
grid
xlabel('Time (ms)');
ylabel('Current (nA)');
legend("Spike-triggered average")

%% Compute a nonlinearity
% (function relating the projection of the injected current onto the STA vs. the firing rate)

% Compute the linear projection of the stimulus onto the STA
fprintf('Computing the nonlinearity...\n')
u = zeros(length(current), 1);
for t=(2*tau):(length(time)-tau)
  if mod(t, 1000/dt) == 0
    fprintf('%f%% done\n', 100 * t / length(time));
  end
  current_slice = current(t-2*tau+1:t+tau)';
  u(t) = mean(current_slice .* sta);
end
%%
% Compute the nonlinearity by estimating E(rate | u)
% bin the spike times using the time array
binned_spikes = hist(spike_times(:, 2), time);

% nonlinearity via a ratio of histograms
nbins = 240;
bins = linspace(-450, 450, nbins);
raw = hist(u, bins);
raw = raw / sum(raw);
conditional = hist(u(spike_times(:,1)), bins);
conditional = conditional / sum(conditional);

% compute the nonlinearity (ratio of histograms)
nln = mean(binned_spikes) * conditional ./ raw;

figure(2);
subplot(2, 2, 2)
plot(bins, raw, 'ro-');  hold on;
set(gca, 'YGrid', 'off', 'XGrid', 'on')
plot(bins, conditional, 'bo-');  hold off;
xlim([-60,60])
ylim([0,0.35])

legend("P(STA projection)", "P(STA projection|spike)", 'location', 'northwest')
xlabel("STA projection")
title("")
subplot(2,2, 4)
minraw = 1e-5;
plot(bins(raw > minraw), nln( raw > minraw) / dt, 'k.-');
xlabel("STA projection")
ylabel("Relative firing rate")
set(gca, 'YGrid', 'off', 'XGrid', 'on')
legend("nonlinearity", 'location', 'northwest')
xlim([-60,60])
ylim([0,0.8])
sgtitle("High-variance input")
