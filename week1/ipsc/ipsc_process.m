datadir = 'C:\Users\DataAnalysis\Vasily\nepr208\1\ipsc\';
threshold = -20;
start_time = 300 + 10;
end_time = 1795;

%%
files = dir([datadir, '*.txt']);

%%
kinetics = zeros(size(files));
npeaks = zeros(size(files));
overshoots = zeros(size(files));
voltagedata = zeros([size(files, 1), 18000]);

for f = 1:size(files)
    fname = files(f).name;
    data = read_EOTN_file([datadir, fname]);

    %
    [py, px] = findpeaks(data.Voltageconcentrationt, data.time);
    %threshold = quantile(data.Voltageconcentrationt, 0.999);

    baseline = median(data.Voltageconcentrationt((end_time < data.time)));
    
    kinetics(f) = str2double(fname(1:(end-4)));    
    npeaks(f) = sum(py > threshold);
    overshoots(f) = max(data.Voltageconcentrationt(start_time < data.time) - baseline);
    voltagedata(f,:) = data.Voltageconcentrationt;
    
    disp(f)
end

%%

subplot(2,3,1)
plot(data.time, voltagedata(10,:))
hold on
plot(data.time, voltagedata(3,:))
plot(data.time, voltagedata(1,:))
xlabel('time, ms')
ylabel('Voltage, mV')
legend(["IPSC kinetics multiplier = 2^{4}", "IPSC kinetics multiplier = 2^{-3}", "IPSC kinetics multiplier = 2^{-5}"])
hold off
xlim([200,1200])


subplot(2,3,4)
plot(data.time, voltagedata(10,:))
hold on
plot(data.time, voltagedata(3,:))
plot(data.time, voltagedata(1,:))
xlabel('time, ms')
ylabel('Voltage, mV')
legend(["IPSC kinetics multiplier = 2^{4}", "IPSC kinetics multiplier = 2^{-3}", "IPSC kinetics multiplier = 2^{-5}"])
hold off
xlim([200,1200])
%%
subplot(1,3,2)
semilogy(log2(kinetics), overshoots, 'o-')
xlim([-5,6])
set(gca,'xtick',-5:1:6)
xtick = get(gca, 'XTick');


str = cellstr( num2str(xtick(:),'2^{%d}') );
format_ticks_v2(gca,str, ' ')

xlabel(["", "IPSC kinetics multiplier"]);
ylabel("Overshoot, mV") 



subplot(1,3,3)
stem(log2(kinetics), npeaks, '')
set(gca,'xtick',-5:1:6)
xlim([-5,6])
xtick = get(gca, 'XTick');
str = cellstr( num2str(xtick(:),'2^{%d}') );
format_ticks_v2(gca,str, ' ')

xlabel(["", "IPSC kinetics multiplier"]);
ylabel("# of action potentials") 

set(gcf,'color','w');
