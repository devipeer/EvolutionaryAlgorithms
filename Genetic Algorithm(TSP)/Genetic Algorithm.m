%% GA main
clear all;
% Parameters:
P=100% 250; %100, 300, 500 % number of chromosomes
n=0.5% 0.8; % 0.5, 0.7, 0.9 percentage of the population that will become parents
pm=0.1% 0.2; % 0.1, 0.7, 0.9 mutations rate
Tmax=1000; % number of iterations
N=10; % number of cities

% Initial population:
P0=[1 2 3 4 5 6 7 8 9 10];
Population=initPop(P0,P);
% Distance matrix for map number 1:
x = [0 3 6 7 15 12 14 9 7 0];
y = [1 4 5 3 0 4 10 6 9 10];
plot(x,y,'black*','MarkerSize',4,'LineWidth',3);
title('Map of cities (number 1)'); hold on;
xlabel('Coord x');
ylabel('Coord y');
xlim([min(x)-2,max(x)+2]);
ylim([min(y)-2,max(y)+2]);
labels = string(1:10);
text(x-.1,y-.3, labels);
% adding connections between different points:
connect=nchoosek(1:length(x),2)';
plot(x(connect),y(connect),"Color",'#FFC0CB');hold on;
% compute the distance matrix [10x10]:
Distance=zeros(length(P0)); %initialize distance matrix
for i=1:N
    for j=1:N
        Distance(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
    end
end
% main itertion:
for i=1:Tmax

    % Evaluation/ cost value:
    [fi,di]=evaluation(Population, Distance);
    mf=max(fi);
    ti=(mf-fi);
    ts=sum(ti);

    % select of parent population: n*P
    PPop=n*P ;% Quantity of parents chromosomes

    for i=1:PPop
       roulette=0;
       random=rand()*ts;
       dd=randperm(P);
       % first individual who exceed the random value is taken to the parents
        for k=1:P
            roulette=roulette+ti(dd(k));
            if roulette>=random
                Parents(i,:)=Population(dd(k),:);
            end
        end
    end
    half=length(Parents)/2;
    P1=Parents(1:half,:);
    P2=Parents(half+1:length(Parents),:);

    %Cycle crossover
    for i=1:length(P1)
        Offspring1(i,:)=cycle_cross(P1(i,:),P2(i,:));
        Offspring2(i,:)=cycle_cross(P2(i,:),P1(i,:));
    end
    Offspring=[Offspring1;Offspring2];

    % Mutation
    for j=1:length(Offspring)
        random_pm= rand();
        if random_pm<pm
             pos_mut=randperm(N,2); % random position of mutations
             temp=Offspring(j,pos_mut(1));
             Offspring(j,pos_mut(1))=Offspring(j,pos_mut(2));
             Offspring(j,pos_mut(2))=temp;
        end
    end
    % Evaluation of Offspring
    [fi_O,di_O]=evaluation(Offspring, Distance);
    % %Join offspring population with original population
    Mixed=[Population;Offspring];

    %Create new population (250 chromosomes) which have min. fi from mixed population
    Mixed_fi=[fi fi_O];
    [values, index]=mink(Mixed_fi,P);
    New_population=Mixed(index,:);
    fi_new=Mixed_fi(index);
    Population=New_population;
end
% Show the best distance and sequence:
[value1, index1]=min(fi_new);
Best_way=Population(index1,:)
Minimal_fi=value1
plot(x([Best_way,Best_way(1)]),y([Best_way,Best_way(1)]),"k.", 'MarkerSize',2,'LineWidth',1)

function Population = initPop(P0,P)
    PP=perms(P0);
    [x2 , ~]=size(PP);
    Pop=[];
    for i=1:P
        Pop1= PP(randi(x2),:);
        Pop=[Pop;Pop1];
    end
    Population=Pop;
end

function [fi,di]=evaluation(Pi,D)
    % Pi = population of P chromosome
    % D= matrix of distance beetween every cities
    [xP, yP]=size(Pi);
    di=zeros(xP,yP);
    for j=1:xP
        tour=0;
        %Distance from the last city till first city:
        di(j,1)=D(Pi(j,1),Pi(j,yP));
        for k=1:(yP-1)
            di(j,k+1)=D(Pi(j,k),Pi(j,k+1)); % distance between every point in given sequence
            tour=tour+di(j,k+1);
        end
        f(j)=di(j,1)+tour; % Total distance in a given sequence
    end
    fi=f;
    di=di;
end
function O=cycle_cross(P1,P2)
    O=[];
    [~, y1]=size(P1);
    O1 =zeros(1,y1);
    O1(1)=P1(1); pos1=1;
    for i=1:length(P1)
        pos1 = find(P1==P2(pos1));
        O1(pos1) = P1(pos1);
        if pos1==1
            empty=find(O1==0);
            O1(empty)=P2(empty);
        end
    end
    O=[O;O1];
    OO=O;
end