clear all
model = readmatrix("model5.txt");

sum_err = 0;
P_size = 100;
lambda = 5;
new_pop = zeros(P_size,7);

method = 2; % 1 - lambda, 2 - lamda + u

%initialize first population
tic
for i=1:P_size
    new_pop(i,1:3)=-10+rand(1,3)*20;
    new_pop(i,4:6)=rand(1,3)*10;
end

%evalute first population
new_pop = evaluate(new_pop, model);
epsilon = 0;
g=0;
while g < 150
    PR = new_pop(1,7);
    mutants = repmat(new_pop, lambda, 1);
    mutants = mutation(mutants);
    mutants = evaluate(mutants,model);
    switch method
        case 1
            new_pop = mutants;
            new_pop = sortrows(new_pop,7);
        case 2
            new_pop = [mutants; new_pop];
            new_pop = sortrows(new_pop,7);
    end
    new_pop = new_pop(1:P_size,:);
    PP = new_pop(1,7);
    epsilon = abs(PR - PP);
    if epsilon < 10^(-5) && epsilon > 0
       break;
    end
    g = g+1;
end
toc

g
figure
plot(model(:,1),model(:,2),".")
hold on 
a = new_pop(1,1);
b = new_pop(1,2);
c = new_pop(1,3);
plot(model(:,1),a*(model(:,1).^2 - b*cos(c*pi*model(:,1))))

function [mutants] = mutation (Pop_to_mut)
    tau1 = 1/sqrt(12);
    tau2 = 1/sqrt(2*sqrt(6));
    m_size = length(Pop_to_mut);
    mutants = zeros(m_size,7);
    for i = 1:m_size
        mutants(i,1) = Pop_to_mut(i,1) + randn()*Pop_to_mut(i,4);
        mutants(i,2) = Pop_to_mut(i,2) + randn()*Pop_to_mut(i,5);
        mutants(i,3) = Pop_to_mut(i,3) + randn()*Pop_to_mut(i,6);
        r_sigma1 = tau1 * normrnd(0,1);
        r_sigma2= tau2 * normrnd(0,1);
        sigmaA = Pop_to_mut(i,4) * exp(r_sigma1) * exp(r_sigma2);
        r_sigma2= tau2 * normrnd(0,1);
        sigmaB = Pop_to_mut(i,5) * exp(r_sigma1) * exp(r_sigma2);
        r_sigma2= tau2 * normrnd(0,1);
        sigmaC = Pop_to_mut(i,6) * exp(r_sigma1) * exp(r_sigma2);
        mutants(i,4:6) = [sigmaA, sigmaB, sigmaC];
    end
end

function [evalueted] = evaluate(Pop_to_ev, model)
    P_size = length(Pop_to_ev);
    evalueted = Pop_to_ev;
    for l=1:P_size
        sum_err = 0;
        a = evalueted(l,1);
        b = evalueted(l,2);
        c = evalueted(l,3);
        for t=1:101
            i = model(t,1);
            o = a*(i.^2 - b*cos(c*pi*i));
            err = (o - model(t,2)).^2;
            sum_err = sum_err + err;
        end
        mid_sq_err = 1/101*sum_err;
        evalueted(l,7) = mid_sq_err;
    end
    evalueted = sortrows(evalueted,7);
end