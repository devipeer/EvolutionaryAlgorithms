x = [3 2 12 7  9  3 16 11 9 2];
y = [1 4 2 4.5 9 1.5 11 8 10 7];
Coords = [x ; y]';
N = size(Coords,1);

%create a distance matrix 
Distances = zeros(N,N);
for i=1:N
    for j=1:N
        Distances(i,j)= ((x(i) - x(j)).^2 + (y(i) - y(j)).^2).^(1/2);
    end
end

alpha = 1;
beta = 5;
p = 0.5;
m = N;
tmax = 100;
dmax = max(max(Distances));
Tau0 = 1/dmax;
Tau = Tau0*(ones(N,N)-diag(ones(1,N)));

for t=1:tmax
    ants = zeros(m,N);
    ants(:,1) = randi([1 N], [m 1]); % start point
    % Order of visting cities
    for i=1:N-1          
        ants(:,i+1) = selectCity(ants(:,1:(i)), m, alpha, beta, Tau, Distances, N);
    end

    % Calculate whole distance for each ant
    ants_route = zeros(m,1);
    for i=1:m
        for j=1:N-1
            ants_route(i) = ants_route(i) + Distances(ants(i,j),ants(i,j+1));
        end
        ants_route(i) = ants_route(i) + Distances(ants(i,N),ants(i,1));
    end

    % Calculate the total quantity of phermone
    dTau = zeros(N,N);
    for i=1:m
        for j=1:N-1
            dTau(ants(i,j),ants(i,j+1)) = dTau(ants(i,j),ants(i,j+1)) + 1/ants_route(i);
        end
        dTau(ants(i,N),ants(i,1)) = dTau(ants(i,N),ants(i,1)) + 1/ants_route(i); 
    end
    dTau = dTau + dTau';
    Tau = (1-p).*Tau + dTau;
end
index = find(ants_route == min(ants_route));
graph_form(ants(index(1),:), ants_route(index(1)), Coords, N)

function city = selectCity(ants, m, alpha, beta, T, D, N)
    cities = randperm(N, N);
    city = zeros(m,1);
    for i=1:m
        pos = ants(i, end);
        remaining = cities(~ismember(cities,ants(i,:)));
        A = T(pos,remaining).^alpha .* 1./D(pos,remaining).^beta; % Desicion table
        A = A./sum(A); % Probability to go from city i to city j
        if length(A) == 1
            city(i) = remaining(1);
        else
            r = rand;
            threshold = 0;
            j = 1;
            while threshold <= r
                threshold = threshold + A(j);
                j = j + 1;
            end
            city(i) = remaining(j-1);
        end
    end
end

function graph_form(route,distance,Coords,N)
%graphic representation of solution

    x = Coords(:,1);
    y = Coords(:,2);
    figure
%towns position with their 'names'
        scatter(x,y,'o r','filled')
        Names = 1:N; 
        Names = cellstr(num2str(Names(:))); 
        text(x, y, Names)
%all paths between cities
    for i=1:(N-1)
        for j=(i+1):N
             line([x(i),x(j)], [y(i),y(j)], 'LineWidth', 0.2, ...
                'LineStyle', ':');
        end
        r = route([i i+1]);
        line([x(r(1)),x(r(2))], [y(r(1)),y(r(2))], 'LineWidth', 0.8,...
            'Color', 'magenta');
    end
    r = route([1 N]);
    line([x(r(1)),x(r(2))], [y(r(1)),y(r(2))], 'LineWidth', 0.8,...
            'Color', 'magenta');
%solution description
    title(['Ant Systems solution' newline ...
            'Map of cities with the best route marked' newline ...
            'The minimal tour distance is equal to ' num2str(distance)]);
    axis([min(x)-1 max(x)+1 min(y)-1 max(y)+1])
%     saveas(gcf,'AS_TSP_solution.png')
end


