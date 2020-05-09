%% Computes the STA for a ganglion cell
%  [[ VASILY KRUZHILIN ]]
%  [[ MAY 9, 2020 ]]


%% Load the data from the hdf5 file (time and stimulus are sampled every 10ms)
dt = 0.01;  % stimulus sampling rate (seconds)
time = h5read('rgc_data.h5', '/time');
stimulus = h5read('rgc_data.h5', '/stimulus');
spike_times = h5read('rgc_data.h5', '/spike_times'); %% spikes can happen faster than 10 ms, binnig leads to lost time presision

% stimulus = stimulus(randperm(size(stimulus, 1)), randperm(size(stimulus, 2)) );

%% Initialize the STE
% compute the dimensions of the filter (spatial and temporal)
spatial_dim = size(stimulus, 1);         % the number of spatial dimensions in the stimulus
filter_length = 40;                      % the number of temporal dimensions in the stimulus

% cut out the first few spikes that occur before the length of the filter (in seconds) has elapsed
spike_times = spike_times(spike_times > filter_length * dt);

% store the indices of the time array corresponding to each spike
% (you'll use this when computing histograms and the nonlinearity of the LN model)
spike_indices = zeros(length(spike_times), 1);

% initialize an array that will store the spike-triggered ensemble (STE)
% it is a matrix with dimensions given by: the number of spikes and the total of dimensions in the filter

ste = zeros(length(spike_times), filter_length*size(stimulus, 1));


%% Collect stimuli that are part of the STE
fprintf('Collecting the spike-triggered ensemble ...\n');
for i=1:length(spike_times)

  % get the nearest index of this spike time
  % (there are multiple ways you can do this)
  % [[ YOUR CODE HERE ]]
  [~, It] = min(abs(spike_times(i) - time));
  spike_indices(i) = It; 
  % select out the stimulus preceeding the spike, and store it in the `ste` array
  % [[ YOUR CODE HERE ]]
  ste(i, :) = reshape(stimulus(:, (It-filter_length ):(It - 1)), ...
      [filter_length*size(stimulus, 1), 1] );
  
  % update the display
  if mod(i, 1000) == 0
    fprintf('%i of %i\n', i, length(spike_times));
  end

end
fprintf('Done!\n');


%% Compute the spike-triggered average (STA)
%% by averaging over the number of spikes in the STE

sta = reshape(mean(ste, 1), [size(stimulus, 1), filter_length]);

% we'll normalize the STA so that it has unit norm (the scale is arbitrary)
sta = sta / norm(sta);

%% Compute the STC and perform PCA

% We're going  he full set of stimuli causing the nueron
% to spike.

% --------Warm up---------
% 
% 1. How many points are there in this cloud?
%
% 2. What is the dimensionality of this space? (Hint: this is
%     stimulus space, so is should be the same dimensionality as
%     the stimulus
%

% --------Performing PCA--------

% Step 1. Collect stimuli
% We've already done this, and they're in STE

% Step 2. Compute the covariance of the STE
%  (you can use matlab's cov function for this)
stc = cov(ste);% [[ YOUR CODE HERE ]]s

% Step 3. Compute the eigenvalues and eigenvectors of 
%   the covariance (you can use matlab's eig function
%   for this)

[eigvecs, eigvals] = eig(stc);% [[ YOUR CODE HERE ]]

% That's basically it! though we'll add

% Step 4. Visulize
%   this is done in the code below.

% -------------------------------


%% Compute the linear projection of the stimulus onto the STA
u = zeros(length(time), 1);      % the variable `u` will store the projection at each time step
for t=filter_length:length(time) % loop over each time point

  % extract the stimulus slice at this time point
  stimulus_slice = stimulus(:, t-filter_length+1:t);

  % store the linear projection (dot product) of the stimulus slice onto the STA
  
  u(t) = mean(sum(sta.*stimulus_slice, [1,2]));% [[ YOUR CODE HERE ]]

end


%% Compute the histogram of the stimulus projection
% bin the spike times using the time array (hint: use the hist command)

spike_counts = hist(spike_times, time);% [[ YOUR CODE HERE ]]

% discretize the values of u, and get the corresponding bin indices
[n, edges, bin_indices] = custom_histcounts(u, 50);
ub = unique(bin_indices);
nonlinearity = zeros(max(ub)-1,1);
bin_centers = edges(1:end-1) + 0.5 * diff(edges);
for bin_index=1:(max(ub)-1)
  % here, you will need to compute the mean spike count (spike_counts) conditioned
  % on times when the filtered stimulus (u) is within a certain bin
  
  nonlinearity(bin_index) = mean(spike_counts(edges(bin_index) < u & u < edges(bin_index + 1))); % [[ YOUR CODE HERE ]]
end


%% Compute the nonlinearity via a ratio of histograms
%  (this is another way to compute the nonlinearity, it
%   is included here for your reference).
bins = linspace(-6, 6, 30);
raw = hist(u, bins);
raw = raw / sum(raw);       % p(stimulus)
conditional = hist(u(spike_indices), bins);
conditional = conditional / sum(conditional);  % p(stimulus|spike)
nln = mean(spike_counts) * conditional ./ raw; % p(spike|stimulus)

%% Visualization
%  **You only need to add titles and axis labels**
figure(1);
load('colormap.mat');
range = [-0.2, 0.2];
imagesc(sta, range);
colormap(cmap);
colorbar;
axis('image');
title('Spike-triggered average')

figure(2);
plot(diag(eigvals), '.');
xlim([0 size(ste, 2)+1]);
xlabel('Direction in space');
ylabel('Spread');
title('Eigenspectrum');
legend('Original', 'Shuffled')

figure(3);
subplot(131);
imagesc(reshape(eigvecs(:,end), spatial_dim, filter_length), range);
colormap(cmap);
axis('image');
title("PC-1");
subplot(132);
imagesc(reshape(eigvecs(:,end-1), spatial_dim, filter_length), range);
colormap(cmap);
axis('image');
title("PC-2");
subplot(133);
imagesc(reshape(eigvecs(:,end-3), spatial_dim, filter_length), range);
colormap(cmap);
axis('image');
title("PC-3");

figure(4);
plot(bin_centers, nonlinearity / dt, 'ro');
hold on
plot(bins, nln / dt, 'k-');
xlabel('STA response')
ylabel('mean spike count')
