clear all;

% Define input file 
datafile="data/pion_correl_all_data_ascii.dat";

% Variables
N = 2224; % Size of dataset
n = 64;   % Amount of time slices 
convFactor = 197.32 / 0.03505  % Gives units in MeV (hc/a)

% Read the correlation function
in = importdata(datafile,' ',1);
t=in.data(1:n, 1);
tfull=in.data(1:2*n, 1);
tfulleff=in.data(1:2*n-1, 1);

c=[];
cfull=[];
for i = 2:2225
    c(:, i-1) = in.data(1:n, i);
    cfull(:, i-1) = in.data(1:2*n, i);
end

% Apply jackknife resampling
y = [];
cm = mean(c, 2);
for i = 1:size(c, 2)
    y(:, i) = (N*cm - c(:, i))/(N-1);
end
yfull = [];
cmfull = mean(cfull, 2);
for i = 1:size(cfull, 2)
    yfull(:, i) = (N*cmfull - cfull(:, i))/(N-1);
end


figure(1)
subplot(2,1,1);

plot(tfull,log(cfull))
xlim([0 127])
title('All datasets')
xlabel('t') 
ylabel('log(c(t))') 

subplot(2,1,2); 
plot(tfull,log(cmfull))
xlim([0 127])
title('Mean of all datasets')
xlabel('t') 
ylabel('log(c(t))') 

sgtitle('Raw data from simulation')

% Compute the log
y_log = log(y);
y_logm = mean(y_log, 2);

% Compute the effective mass
E_0 = log( y(1:63, :) ./ y(2:64, :) );
E_0full = log( yfull(1:127, :) ./ yfull(2:128, :) );
E_0m = mean(E_0, 2) ;
E_0mfull = mean(E_0full, 2) ;

figure(2)
plot(tfulleff,E_0mfull)
xlim([0 126])
title('Effective mass')
xlabel('t') 
ylabel('m_{eff}') 



% Computing line constants for 10 different intervals 11:36 - 20:45
res = []
for i = 1:size(c, 2)
    for l = 1:10
        interval = 5+l:30+l ;
        [b, a, chi] = linRegress(t(interval), y_log(interval, i));
        res(:, l, i) = [b, a, chi];
    end
end

% Getting the mean of the results
res_mean = mean(res, 3);

% Calculating statistical error by std
err_stat = std(mean(res(1, :, :), 2)) * convFactor
err_stat2 = sqrt((N-1)/N ) 
% Calculating systematic error by half observation difference
err_sys = 0.5 * (max(res_mean(1, :)) - min(res_mean(1, :))) * convFactor
% Calculation effective mass
E = mean(res_mean(1, :)) * -convFactor

% Calculation errors for error bars 
err_E = sqrt( (N-1)/N * sum( (E_0 - E_0m * ones(1, 2224)).^2, 2) );

% Plots
figure(3)
plot(t, y_logm)
figure(4) 
plot(res_mean(1, :));
title('b')
figure(5) 
plot(res_mean(2, :));
title('a')
figure(6) 
plot(res_mean(3, :));
title('chi^2')

figure(7)
errorbar(t(1:63), E_0m, err_E) 

