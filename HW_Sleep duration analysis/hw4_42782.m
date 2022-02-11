% % Q1 a)
load('sleep.mat')
x = (1:1:180);
scatter(x,sleephours, 24);
xlim([1, 180])
xlabel('Time(day)')
ylabel('Sleep Duration(hr)')

% % Q1 b)
ave_sleep = roundn(sum(sleephours)/180, -2);
c = xcorr(sleephours-ave_sleep,'coeff');
plot(x, c(180:end,1));
xlim([1, 180])
xlabel('Time(day)')
ylabel('Sleep Duration Autocorrelation')

% Q1 e) sleep exactly 7 hours
t_fix = ones(180,1)*7;
plot(x, t_fix);
[c_fix, lag_fix] = xcorr(t_fix-7,'coeff');
plot(x, c_fix(180:end,1));
xlim([1, 180]);
ylim([-1, 1]);
legend('Correlation');
xlabel('Time(day)')
ylabel('Fixed Sleep Duration Autocorrelation')


% % Q2 a)
load('pain.mat')
sec = (0:0.5:699.5);
[AX,H1,H2]=plotyy([sec',sec'],[pu_left,pu_right],sec',therm);
set(AX,'Xlim',[0,700])
set(AX(1),'Ylim',[0, 2]) 
set(AX(2),'Ylim',[20, 60]) 
set(get(AX(1),'Ylabel'),'string',{'peripheral blood flow from the left forearm)';'peripheral blood flow from the right forearm'});
set(get(AX(2),'Ylabel'),'string','heat pain pulses (°C)');
xlabel({'Time (sec)'});
legend([H1(1),H1(2),H2],'Left','Right','pulse');

% Q2 b) i, ii
[c_pur, lag_pur]= xcorr(therm-(sum(therm)/1400), pu_right-(sum(pu_right)/1400), 'coeff');
time = lag_pur/fs;
plot(time, c_pur);
xlabel('Time (sec)');
ylabel('Cross-correlation');
title('CC between heat pain pulses and right forearm-peripheral blood flow');
[pri, mri] = min(c_pur);
[prx, mrx] = max(c_pur);

[c_pul, lag_pul] = xcorr(therm-(sum(therm)/1400), pu_left-(sum(pu_left)/1400), 'coeff');
time = lag_pul/fs;
plot(time, c_pul);
ylim([-0.65 0.3]);
xlabel('Time (sec)');
ylabel('Cross-correlation');
title('CC between heat pain pulses and left forearm-peripheral blood flow');
[pli, mli] = min(c_pul);
[plx, mlx] = max(c_pul);


