%% 
% Number of individuals
N = 100; %total number of individuals in the population
NH1 = 99; %healthy w/o immunity on day 1
NI1 = 1; %infected on day 1
NC1 = 0; %contagious on day 1
NHI1 = N - (NI1 + NC1 + NH1); %healthy with immunity on day 1

% Model parameters
Tinc = [3, 14]; %incubation period [days]
Trec = [7, 50]; %recuperation period [days]
Timm = [2*365, 10*365]; %immunity period [days]
Nsoc = 0.2*N; %average number of times of people making contact per day
SDsoc = 0.1*Nsoc; %standard deviation of the number of times of people making contact per day
Pbio = 0.3; %probability of contracting disease if contact occurs
nday = 360; %number of simulation days

NH = zeros(nday,1); %no. of indiv who are healthy w/o immunity on each day
NI = zeros(nday,1); %no. of indiv who are infected on each day
NC = zeros(nday,1); %no. of indiv who are contagious on each day
NHI = zeros(nday,1); %no. of indiv who are healthy w/ immunity on each day

NH(1) = NH1;
NI(1) = NI1;
NC(1) = NC1;
NHI(1) = NHI1;

HS = [zeros(NH1,1);... %no. of healthy w/o immunity
-randi(Tinc, NI1, 1);... %no. of infected
randi(Trec, NC1, 1);... %no. of contagious
zeros(NHI1, 1)];%no. of healthy w immunity

IM = [zeros(NH1,1);...
zeros(NI1,1);...
zeros(NC1,1);...
randi(Timm, NHI1, 1)];

Nexposure = max(0,floor(normrnd(Nsoc, SDsoc))); % exposures/traffics per day (integer)
time = linspace(1, nday, nday);
nsim = 100;
NH_all = zeros(nday,nsim);
NI_all = zeros(nday,nsim);
NC_all = zeros(nday,nsim);
NHI_all = zeros(nday,nsim);
for g = 1: nsim
    [NH,NI,NC,NHI] = spreadDisease(HS, IM, nday, NH1, NI1, NC1, NHI1, Tinc, Trec, Timm, Nexposure, Pbio, N);
    NH_all(:,g) = NH;
    NI_all(:,g) = NI;
    NC_all(:,g) = NC;
    NHI_all(:,g) = NHI;
end
figure;
plot(time, NH_all, 'b',  time, NI_all, 'g', time, NC_all, 'r', time, NHI_all, 'm');
handle = get(gca,'Children');
legend([handle(3*nsim+1), ...
handle(2*nsim+1), ...
handle(nsim+1), ...
handle(1)],...
'Healthy', ... 
'Infected', ... 
'Contagious', ... 
'Healthy IM');


subpop_all = {NH_all(180, :);  NI_all(180, :); 
           NC_all(180, :); NHI_all(180, :)}; 
subpop_label = {'Healthy w/o imm';  'Infected';
             'Contagious';    'Healthy w/i imm'}; 
subpop_color = {'b', 'g', 'r', 'm'};
npop = length(subpop_all); 
figure; 
for k=1:npop 
    ax(k)=subplot(1,4,k); histogram(subpop_all{k},'FaceColor', subpop_color{k}); 
    title(subpop_label{k}); 
    ylabel('# of Simulation'); 
    xlabel('Cases on day 180'); 
end 

function [NH,NI,NC,NHI] = spreadDisease(HS, IM, nday, NH1, NI1, NC1, NHI1, Tinc, Trec, Timm, Nexposure, Pbio, N)
    NH = zeros(nday,1);
    NI = zeros(nday,1); 
    NC = zeros(nday,1); 
    NHI = zeros(nday,1); 
    
    NH(1) = NH1;
    NI(1) = NI1;
    NC(1) = NC1;
    NHI(1) = NHI1;
    for j = 2: nday
        for e = 1 : Nexposure
    
            r = randperm(N, 2);
            i = r(1);
            k = r(2);
            DT = floor(rand + Pbio);
            if DT == 1
                if HS(i)>0 && HS(k)==0 && IM(k)==0
                    HS(k)=-randi(Tinc, 1, 1);
                end
                if HS(i)==0 && IM(i)==0 && HS(k)>0
                    HS(i)=-randi(Tinc, 1, 1);
                end
            end
        end
        NH(j) = sum(HS == 0 & IM == 0);
        NI(j) = sum(HS < 0 & IM == 0);
        NC(j) = sum(HS > 0 & IM == 0);
        NHI(j) = sum(HS == 0 & IM < 0);
        for n=1:N
            if IM(n)<0  
                IM(n) = IM(n) + 1;
            elseif HS(n) < 0 && HS(n) ~= -1
                HS(n) = HS(n) + 1;
            elseif HS(n) == -1
                HS(n)=randi(Trec, 1, 1);
            elseif HS(n) > 0 && HS(n) ~= 1
                HS(n) = HS(n) - 1;
            elseif HS(n) == 1
                HS(n) = 0;
                IM(n)=-randi(Timm, 1, 1);
            end
        end
    end
end