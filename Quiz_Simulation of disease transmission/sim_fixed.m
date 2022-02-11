N = 10; %total number of individuals in the population
NH1 = 9; %healthy w/o immunity on day 1
NI1 = 1; %infected on day 1
NC1 = 0; %contagious on day 1
NHI1 = N - (NI1 + NC1 + NH1); %healthy with immunity on day 1

% Model parameters
Tinc = [3, 3]; %incubation period [days]
Trec = [7, 7]; %recuperation period [days]
Timm = [100, 100]; %immunity period [days]
Nsoc = N; %average number of times of people making contact per day
SDsoc = 0; %standard deviation of the number of times of people making contact per day
Pbio = 1; %probability of contracting disease if contact occurs
nday = 25; %number of simulation days

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
i_all = [ 4; 1; 7;  4; 4;  6; 8; 2; 5; 1];   
k_all = [10; 2; 5; 10; 3; 10; 6; 6; 7; 2];
for j = 2: nday
    for e = 1 : Nexposure
        i = i_all(e);
        k = k_all(e);
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

% Plot
time = linspace(1, nday, nday);
figure;
plot(time,NH, time,NI, time,NC, time,NHI);
xlabel('Time (days)');
ylabel('Number of individuals in subpopulation');
legend('NH', 'NI', 'NC', 'NHI');


function [NH,NI,NC,NHI] = spreadDisease(HS, IM, nday, NH1, NI1, NC1, NHI1, Tinc, Trec, Timm, Nexposure, Pbio)
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
                if HS(i)==0 && IM(i)==0 && HS(k)>0 
                    HS(i)=-randi(Tinc, 1, 1);
                end
                end
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
        elseif HS(n) < -1
            HS(n) = HS(n) + 1;
        elseif HS(n) == -1
            HS(n)=randi(Trec, 1, 1);
        elseif HS(n) > 1
            HS(n) = HS(n) - 1;
        elseif HS(n) == 1
            HS(n) = 0;
            IM(n)=-randi(Timm, 1, 1);
        end
    end
end
