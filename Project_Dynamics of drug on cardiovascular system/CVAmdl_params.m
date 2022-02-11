% stroke volume/L
SV = 0.097;  
% Right atrial pressure/mmHg
RAP = 5;
% BR sensitivty to change in MAP
S_p = 0.01;
% BR sensitivty to the rate of change in MAP
S_r = 0.01;
% Relative HR baroreflex gain 
alpha = 1;
% Relative TPR baroreflex gain
beta = 1;
% Time constant ~ total baroreflex dynamic/hr
t = 0.01;
% Time constant ~ HR baroreflex dynamic/hr
t1 = 0.01;
% Time constant ~ TPR baroreflex dynamic/hr
t2 = 0.01;
% Pre-infusion equilibrium of MAP/mmHg
MAP_eq = 90;
% Pre-infusion equilibrium of HR/bpm
HR_eq = 70;
% Pre-infusion equilibrium of TPR/RU+
TPR_eq = 12.5;
% Drug infusion period/hr
T = 10;
% Plasma volume of drug distribution/L
V = 10;
% Elimination rate constant/hr-1
gamma = 0.7;
% Drug infusion rate/μg∙hr-1
R0 = 350;
% Pharmacodynamic gain factor/RU∙L∙μg-1
m = -0.09;
