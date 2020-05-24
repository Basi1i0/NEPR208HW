% Explorations with the Hopfield Model
% Surya Ganguli June 3rd 2014

close all;

N = 2000; % number of neurons
P = 100;  % number of patterns to be stored
T = 300;  % maximial number of iterations to run Hopfield Model

% Note when N=1000, the Hopfield capacity is P ~ 0.14N ~ 140.  Thus the
% network should not be able to recall stored memories when P > ~ 140.



K = 20;  % number of spins/neurons to flip in scenario InitCondType = 2;

Kmax = N*0.7;
dK = N/100;
Ks = (0:dK:Kmax) + 1;

Pmax = N*0.6;
dP = N/100;
Ps = (0:dP:Pmax) + 1;

nK = length(Ks);
nP = length(Ps);
nrep = 120;

%%
overlaps = nan(nK, nP, nrep);

parfor rep = 1:nrep
    disp(rep/nrep);
    
    for iP = 1:nP
        P = Ps(iP);
        patterns = sign(randn(N,P));  % patterns(:,i) = column vector of the i'th random pattern
        J = sqrt(1/N) * patterns * patterns';  % Hopfield connectivity matrix to store patterns
        InitCondType = 2;  % type of initial condition
                              % 1 = completely random
                              % 2 = K random flips from a randomly chosen stored pattern
        for iK = 1:nK
            K = Ks(iK);

            if InitCondType == 1
                InitCond = sign(randn(N,1));
            else
                PattChoice = ceil(P*rand(1,1));
                InitCond = patterns(:,PattChoice);
                NeuronChoice = [ones(K,1) ; zeros(N-K,1)];
                NeuronChoice = find(NeuronChoice(randperm(N))>0);
                % Now NeuronChoice is a vector of K unique random numbers from 1 to N 
                % corresponding to K neurons to flip
                InitCond(NeuronChoice) = -InitCond(NeuronChoice);  % flip the K randomly chosen neurons
            end

            % Run the Hopfield Dynamics and store the network trajectory
            Trajectory = zeros(N,T);  
            Trajectory(:,1) = InitCond;
            iter = 2;
            scurr = InitCond;
            sprev = zeros(N,1);
            while (sum((scurr-sprev).^2) > 0) & (iter <= T)  % keep going unless you have hit a fixed point or reached T iterations
                sprev = scurr;
                scurr = sign(J*scurr ); %+ 10*randn(length(J*scurr), 1) 
                Trajectory(:,iter) = scurr;
                iter = iter + 1;
            end

            NumIter = iter-1;  % This is now either the number of iterations it took to get to
                             % a fixed point, or T if it never reached a fixed point
        % The overlap between two
        % patterns sa and sb is 1/N \sum_i sa_i * sb_i .  It is always a number between
        % -1 and 1 with 1 indicating perfect overlap
        %
        % A pattern is recalled if the final fixed point has a large overlap with
        % that pattern (say absolute value of overlap greater than 0.9) and small overlap with the rest
            overlap = (1/N) * patterns' * Trajectory;

            % determine whether or not recall was successful
            successful_recall = max(overlap(:,NumIter)) > 0.99;

            overlaps(iK, iP, rep) = max(overlap(:,NumIter)) ;
            end
    end
end
%%
threshold = 0.9;
subplot(1,3,1)
imagesc(Ks/N, Ps/N, mean( abs(overlaps1) > threshold, 3));
ylabel("proportion of initial flips K/N");
xlabel("relative number of stored patterns P/N");
title("Recall probability (threshold = 0.9)")
colormap(mm.redblue)
colorbar

threshold = 0.99;
subplot(1,3,2)
imagesc(Ks/N, Ps/N, mean( abs(overlaps1) > threshold, 3));
ylabel("proportion of initial flips K/N");
xlabel("relative number of stored patterns P/N");
title("Recall probability (threshold = 0.99)")
colormap(mm.redblue)
colorbar

threshold = 0.999;
subplot(1,3,3)
imagesc(Ks/N, Ps/N, mean( abs(overlaps1) > threshold, 3));
ylabel("proportion of initial flips K/N");
xlabel("relative number of stored patterns P/N");
title("Recall probability (threshold = 0.999)")
colormap(mm.redblue)
colorbar


%%
threshold = 0.999;
subplot(1,3,1)
imagesc(Ks/N, Ps/N, mean( abs(overlaps1) > threshold, 3));
ylabel("proportion of initial flips K/N");
xlabel("relative number of stored patterns P/N");
title("Recall probability (N = 1000)")
colormap(mm.redblue)
colorbar

subplot(1,3,2)
imagesc(Ks/N, Ps/N, mean( abs(overlaps2) > threshold, 3));
ylabel("proportion of initial flips K/N");
xlabel("relative number of stored patterns P/N");
title("Recall probability (N = 2000)")
colormap(mm.redblue)
colorbar


subplot(1,3,3)
[~,Ik] = min(abs(Ks-100));
[~,Ip] = min(abs(Ps-200));
mean( abs(overlaps(Ik, Ip, :)) > 0.9, 3)


plot(Ps/N, mean( abs(overlaps1(1, :, :)) > threshold, 3),'o-' ); hold on;
plot(Ps/N, mean( abs(overlaps2(1, :, :)) > threshold, 3),'o-' ); hold off;
xlabel('relative of stored patterns P/N')
ylabel('probability of recall')
title("Recall probability for K = 1")
legend('N = 1000', 'N = 2000')
%%
% Visualize overlap dynamics:  as the network evolves from its initial
% condition, the following will show how the overlap between the trajectory
% and each stored memory evolves over time.  
%
%% UNCOMMENT TO OBSERVE THE NETWORK DYNAMICS
%  (press a key to iterate through the overlap dynamics)
%  (or comment this out to suppress the visualization)
figure
for iter = 1:NumIter
    bar(overlap(:,iter));
    xlim([0 P]);
    xlabel('Pattern Index');
    ylabel('Overlap');
    if iter == NumIter
      title(['Iteration = ', num2str(iter)]);
    else
      title(['Iteration = ', num2str(iter), ', (press any key to continue)']);
    end
    set(gca,'Ylim', [-1,1]);
    drawnow;
    waitforbuttonpress;
end
