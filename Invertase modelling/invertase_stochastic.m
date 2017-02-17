% Stochactic modelling of invertase kinetics
% (c) 2017 by Group 3

% The script performs a Gillespie type stochastic simulation of dimeric
% invertase kinetics, applicable to e.g. PhiC31, Bxb1 and Tp901

% Modelling based on the following research:
% Pokhilko et al. (2016) The mechanism of PhiC31 integrase directionality:
% experimental analysis and computational modelling.
% in Nucleic Acid Research, Vol. 44, No. 15.

% Kinetic parameters based on:
% http://2016.igem.org/Team:ETH_Zurich/Parameters

function invertase_stochastic
% Rate constants - time in minutes, concentration in nmol/dm3
param.Pon      = 1.0;     % Promoter activity
param.kmRNAint = 0.3465;  % mRNA translation rate (assuming value from practicals)
param.kBxb1    = 7.921;   % Bxb1 transcription rate (assuming value from practicals)
param.kDBxb1   = 1.0; 	  % DBxb1 dimerization rate
param.k_DBxb1  = 10.0;    % DBxb1 dissociation rate
param.kattBP   = 70.0;    % Affinity of DBxb1 to attB and attP binding
param.k_attBP  = 10.0;    % attB/P - DBxb1 dissociation
param.kattLR   = 15.0;    % Affinity of DBxb1 to attL and attR binding
param.k_attLR  = 10.0;    % attL/R - DBxb1 dissociation
param.kflip    = 0.04;    % Switch flipping rate
param.dmRNAint = 0.01995; % mRNA degradation (assuming typical value)
param.dBxb1    = 0.0198;  % Bxb1 degradation (assuming no active degradation pathway)
param.dDBxb1   = 0.0198;  % DBxb1 degradation (assuming no active degradation pathway)

% Initial state
%       Bxb1       S0      S2     attLR0         
%   mRNA    DBxb1      S1      Psw    attLR1
s0 = [0   0   0   100   0   0   0   0   0];      % S0: plasmid copy number

% Transitions
%             Bxb1     S0      S2     attLR0
%         mRNA    DBxb1    S1      Psw    attLR1
process = [ 1   0   0   0   0   0   0   0   0;   % process 01
            0   1   0   0   0   0   0   0   0;   % process 02
            0  -2   1   0   0   0   0   0   0;   % process 03
            0   2  -1   0   0   0   0   0   0;   % process 04
            0   0  -1  -1   1   0   0   0   0;   % process 05
            0   0   1   1  -1   0   0   0   0;   % process 06
            0   0  -1   0  -1   1   0   0   0;   % process 07
            0   0   1   0   1  -1   0   0   0;   % process 08
            0   0   0   0   0  -1   1   0   0;   % process 09
            0   0  -1   0   0   0   0  -1   1;   % process 10
            0   0   1   0   0   0   0   1  -1;   % process 11
           -1   0   0   0   0   0   0   0   0;   % process 12
            0  -1   0   0   0   0   0   0   0;   % process 13
            0   0  -1   0   0   0   0   0   0];  % process 14

% Steps
n = 1000000;

% Random number generator seed
rng(2);

% Run a Gillespie simulation
    [s,t]=gillespie(s0,param,process,n);
    figure;
    plot(t,s(:,7))
    xlabel('Time (minutes)'); ylabel('% Flipped');
    %csvwrite('/PATH/TO/SOME/FOLDER/phic31.csv',[t,s(:,7)]);
end

% Converts current state into process rates
function [rate] = state_to_rate(param,state)
mRNAint   = state(1,1);
Bxb1      = state(1,2);
DBxb1     = state(1,3);
S0        = state(1,4);
S1        = state(1,5);
S2        = state(1,6);
attLR0    = state(1,8);
attLR1    = state(1,9);

rate = [param.kmRNAint*param.Pon;    % process 1
        param.kBxb1*mRNAint;         % process 2
        param.kDBxb1*Bxb1*(Bxb1-1);  % process 3
        param.k_DBxb1*DBxb1;         % process 4
        2*param.kattBP*S0*DBxb1;     % process 5
        param.k_attBP*S1;            % process 6
        param.kattBP*S1*DBxb1;       % process 7
        2*param.k_attBP*S2;          % process 8
        param.kflip*S2;              % process 9
        param.kattLR*attLR0*DBxb1;   % process 10
        param.k_attLR*attLR1;        % process 11
        param.dmRNAint*mRNAint;      % process 12
        param.dBxb1*Bxb1;            % process 13
        param.dDBxb1*DBxb1];         % process 14
end

% Advances the Gillespie simulation by one timestep
function [next_state,dt] = step(state,rate,transition)
% The rate of any type of event occuring
r_total = sum(rate,1);
% Chooses the next event type with uniform distribution
% weighed by the event probabilities rate(i)/r_total
p = random('Uniform',0,1);
count = 0;
    for i=1:size(rate,1)
      p = p - rate(i,1)/r_total;
      if p <= 0
          count = i;
          break;
      end
    end
% Performs the chosen process
next_state = state + transition(count,:);
% Exponential distribution with time constant 1/r_total
dt = random('Exponential',1/r_total);
end

% Runs a Gillespie simulation with given initial parameters
% for n steps. There are no absorbing states due to constant kr.
function [s,t] = gillespie(s0,param,transition,n)
% Initializes matrices to store the evolving time and state
t=zeros(n+1,1);
s=zeros(n+1,size(s0,2));
s(1,:) = s0;
    % Iterate step using most recent state
    for i=1:n
        [s_next,dt] = step(s(i,:),state_to_rate(param,s(i,:)),transition);
        s(i+1,:) = s_next;
        t(i+1,1) = t(i,1) + dt;
    end
end
