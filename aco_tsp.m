% Sample program- Solving Travelling Salesman Problem using ACO
clc; clear;
userview = memory
t=cputime;
% INPUT parameters
% Co-ordinates of towns considered
C       = [0 0; 0.7 0.5; 1 0; 1 1; 0 1];
%m       =5;     % Number of ants
m=size(C,1);
Nc_max	= 50;   % Maximum number of iterations
alpha	= 1;	% Parameter representing the importance of trail
beta	= 5;	% Parameter representing the importance of visibility
rho		= 0.5;	% Evaporation
Q		= 50;	% A constant
n		= size(C,1); % Number of towns
D		= ones(n,n);	% Initializing the Distance array
% Calculation of the distances dij
for i=1:n
    for j=1:n
        if i<j
            D(i,j)	= sqrt((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2);
        end
    D(j,i)			= D(i,j);
	end
end

eta				= 1./D;        % Visibility -which says that close towns should be chosen with high probability
pheromone		= ones(n,n);	% Initializing the pheromeone array
tabu_list		= zeros(m,n);	% List of towns already visited (barred from visiting this town till next iteration)
Nc				= 0;			% Beginning of iteration
routh_best		= zeros(Nc_max,n);
length_best		= ones(Nc_max,1);
length_average	= ones(Nc_max,1);
% Start of iterations
while Nc<Nc_max
    
    rand_position		= [];
    for i=1:ceil(m/n)                   %
        rand_position	= [rand_position,randperm(n)]; % randperm(n) returns a random permutation of the integers 1:n.  
    end
    tabu_list(:,1)		= (rand_position(1:m))';
    for j=2:n
        for i=1:m
            city_visited	= tabu_list(i,1:(j-1));
            city_remained	= zeros(1,(n-j+1));
            probability		= city_remained;
            cr				= 1;
            for k=1:n
                if length(find(city_visited==k))==0
                    city_remained(cr)	= k;
                    cr					= cr+1;
                end
            end
            for k=1:length(city_remained)
                probability(k) 	= (pheromone(city_visited(end),city_remained(k)))^alpha*(eta(city_visited(end),city_remained(k)))^beta;
            end
            probability			= probability/sum(probability);
            pcum				= cumsum(probability);
            select				= find(pcum>= rand);
            to_visit			= city_remained(select(1));
            tabu_list(i,j)		= to_visit;
        end
    end
    if Nc>0
        tabu_list(1,:)			= routh_best(Nc,:); 
    end
   
    total_length				= zeros(m,1);
    for i=1:m
        r						= tabu_list(i,:);
        for j=1:(n-1)
            total_length(i)		= total_length(i)+D(r(j),r(j+1));
        end
        total_length(i)			= total_length(i)+D(r(1),r(n));
    end
    length_best(Nc+1)			= min(total_length);
    position					= find(total_length==length_best(Nc+1));
    routh_best(Nc+1,:)			= tabu_list(position(1),:);
    length_average(Nc+1)		= mean(total_length);
    Nc							= Nc+1;
    delta_pheromone				= zeros(n,n);
    for i=1:m
        for j=1:(n-1)
            delta_pheromone(tabu_list(i,j),tabu_list(i,j+1))	= delta_pheromone(tabu_list(i,j),tabu_list(i,j+1))+Q/total_length(i);
        end
        delta_pheromone(tabu_list(i,n),tabu_list(i,1))			= delta_pheromone(tabu_list(i,n),tabu_list(i,1))+Q/total_length(i);
    end
    pheromone					= (1-rho).*pheromone+delta_pheromone;
    tabu_list					= zeros(m,n);
end
time=cputime-t
position		= find(length_best==min(length_best));
shortest_path	= routh_best(position(1),:)
shortest_length	= length_best(position(1))

figure(1)
set(gcf,'Name','Ant Colony Optimization！！Figure of shortest_path')
N=length(shortest_path);
scatter(C(:,1),C(:,2),50,'filled');
hold on
plot([C(shortest_path(1),1),C(shortest_path(N),1)],[C(shortest_path(1),2),C(shortest_path(N),2)],'r','LineWidth',2)
hold on
xlabel('x-coordinates of cities')
ylabel('y-coordinates of cities')
%axis([-0.25 1.25 -0.25 1.25])
for i=2:N
    plot([C(shortest_path(i-1),1),C(shortest_path(i),1)],[C(shortest_path(i-1),2),C(shortest_path(i),2)],'r','LineWidth',2)
    hold on
end

figure(2)
set(gcf,'Name','Ant Colony Optimization！！Figure of length_best and length_average')
plot(length_best,'ro')
hold on
plot(length_average,'k','LineWidth',2)
xlabel('No. of iterations')
ylabel('Distance')
legend('Minimum distance','Average distance travelled');
memory