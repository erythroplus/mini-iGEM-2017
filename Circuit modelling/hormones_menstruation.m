minute = 1:1:30;
estradiol = [118.42,133.01,133.17,125.95,134.65,151.33,159.48,170.34,195.9,228.2,269.07,343.68,477.27,661.19,914.84,780.76,320.59,261.32,338.94,454.07,499.49,497.07,531.14,504.39,499.18,526.68,350.65,322.24,229.7,249.28];
progesterone = [1.66,1.27,1.02,0.82,0.74,0.94,0.72,0.67,0.73,0.51,0.57,0.59,0.51,0.59,1.02,2.66,5.02,12.13,20.84,29.75,36.08,36.52,40.32,39.65,33.76,34.11,18.24,16.9,10.73,9.5];

% estradiol is in pmol / L 
% progesterone is in nMol / L, if we write it in pico we have more parts
progesterone = progesterone/1000;

%% Progesteron sense 1 design
% Use 40 days 24*60*40
tspan = [0 57600];
[t,y] = ode113(@react_combined,tspan,[10 10 1 1 1 0 10 10 1 1 1 0 0 0 0 0 0 0]);
y(:,5) = (y(:,5) > 0) .* y(:,5);
y(:,11) = (y(:,11) > 0) .* y(:,11);
% Gaussian estimates
f = 7897*exp(-((t-2.746*10^4)/2597).^2) + 3.904*10^4*exp(-((t-3.252*10^4)/6783).^2);
g = 0.7007.*exp(-((t-1.929.*10.^4)./4768).^2) + 0.5464.*exp(-((t-3.294.*10^4)./8017).^2);
subplot(3,1,1);
% Fast conversion to days
t = t./24./60;
p = plotyy(t,y(:,5),t,f);
set(p,'linewidth',2.0); set(gca,'FontSize',14);
grid on;
xlabel('day');
ylabel('concentration (nMol)');
title('Concentration response of TEV protein to progesterone'); 
legend('TEV protein', 'Progesterone hormone');
subplot(3,1,2);
q = plotyy(t,y(:,11),t,g);
set(q,'linewidth',2.0); set(gca,'FontSize',14);
grid on;
xlabel('day');
ylabel('concentration (nMol)');
title('Concentration response of TALA protein to estradiol'); 
legend('TALA protein', 'Estradiol hormone');
subplot(3,1,3);
r = plot(t,y(:,18));
set(r,'linewidth',2.0); set(gca,'FontSize',14);
grid on;
xlabel('day');
ylabel('concentration (nMol)');
title('Sensor response to combined hormone cycle'); 

% All differential equations on a plot
figure;
for i=1:1:18
    subplot(4,5,i)
    plot(t,y(:,i));
end

day = 1:1440:41761;
minute_interval = 1:1:41761;
estradiol = interp1(day,estradiol,minute_interval);
progesterone = interp1(day,progesterone,minute_interval);
figure;
plot(minute_interval, estradiol)
xlabel('minutes');
ylabel('Concentration [nMol]');
title('Estrogen');
hold on;
plot(t, g)
xlabel('minutes');
ylabel('Concentration [nMol]');
legend('Original Estrogen Data (Mean)','Gaussian Curve Fit');

figure;
plot(minute_interval, progesterone)
xlabel('minutes');
ylabel('Concentration [nMol]');
title('Progesteron ');
hold on;
plot(t, f)
xlabel('minutes');
ylabel('Concentration [nMol]');
legend('Original Progesteron Data (Mean)', 'Gaussian Curve Fit');
