randSeed = 11;
inputNeuronNum = 2;
patternNum = 100;
learningRate = 0.1;
epochNumMax = 1000;


rng('default');
rng(randSeed);
%%
f_and = @(s1,s2) ((s1 == 1).*( s2 == 1)*2) -1 ;
%%f_and2 = @(s1,s2) (s1.*s2);
f_xor = @(s1,s2) xor((s1+1)/2,(s2+1)/2)*2 - 1 ;

%%  Define network and connectivity
perceptron = MyPerceptron(inputNeuronNum);

%%  Define patterns
s = sign(randn(inputNeuronNum,patternNum));
sigma = f_xor(s(1,:), s(2,:)); %sign(randn(patternNum,1));

%%
perceptron.Train(s, sigma, learningRate, 1000, true)
%%
scatter(s(1,sigma == 1), s(2,sigma == 1), 'og','filled')
hold on
scatter(s(1,sigma == -1), s(2,sigma == -1), 'or','filled')

set(gca,'xtick',-1:1:1)
set(gca,'ytick',-1:1:1)

xs = -2:0.1:2;
ys = -(xs*perceptron.connVec(2) + perceptron.connVec(3)) ./ perceptron.connVec(1);
plot(xs, ys, '--', 'color', 'black')
legend(' = +1', ' = -1', 'separation boundary')

xlim([-1.1, 1.1])
ylim([-1.1, 1.1])
grid()

%%
ntrials = 1000;
maxNeuronNum = 15;
maxpatterns = zeros(ntrials, maxNeuronNum);

for inputNeuronNum = 1:maxNeuronNum
    for trial = 1:ntrials
        perceptron = MyPerceptron(inputNeuronNum);
        s = permn([-1,1], inputNeuronNum)';
        sigma = sign(randn(size(s,2),1))';
        
        for patternNum = 1:size(s, 2)
            success = ...
                perceptron.Train(s(:, randperm(size(s, 2), patternNum)), sigma(:, 1:patternNum), learningRate, 1000);
            if(~success)
                patternNum = patternNum -1;
                break; 
            end
        end

        maxpatterns(trial, inputNeuronNum) = patternNum;
    end
    disp(inputNeuronNum);
end


%%
mean_values = mean(maxpatterns);
x = 1:maxNeuronNum; 
p = polyfit(x(5:maxNeuronNum), mean_values(5:maxNeuronNum),1);
errorbar(mean(maxpatterns), std(maxpatterns,[],1)/sqrt(ntrials), 'ob', 'MarkerFaceColor', 'b')
hold on
plot(x, p(1)*x + p(2), 'r')
ylabel('Average # of patterns')
xlabel('# of input neurons')
legend('estimated', 'y \approx 1.7 x')
grid()
