%%%%%%%


%%%%%% initialise model parameters

%%%%%%%%%%
N = 120; % size of the domain for finite difference

Ly = 120; % \nu m height of the migratory domain

chemo = ones(N,N); % initial chemoattractant distribution
chemo_new = ones(N,N); % store new concentration of chemical

dt = 1; % time step most likely wrong

D = 0.1; % diffusion coefficient of chemoattractat

t = 0; % initialise time


di = 1; %% space step in x direction WRONG
dj = 1; %% space step in y direction WRONG

kai = 0.0001; %1/h production rate of chemoattractant 

% parameter for internalisation
R = 7.5; % \nu m cell radius
lam = 100 ; % to 1000 /h chemoattractant internalisation rate

%%%%%%%%%% initialise first cells
initial_number_cells = 10;

cells = ones(initial_number_cells,2); % two columns correspond to x and y coordinates

t_final = 100;

%%%%%% angle choice parameters
alpha = 0;
beta = 1;


for i = 1:initial_number_cells
    cells(i,2) = N*rand;
end

%%% generate a figure with initial cells
figure
scatter(cells(:,1),cells(:,2),'b','filled')

intern_term = zeros(N,N); % internalisation term initially zero

%%%%%%% solve chemoattractant profile
while t < t_final
    
    for i = 2:N-1

        for j= 2:N-1

            
            %%%%%% internalisation
            for k = 1: size(cells,1)
                intern_term(i,j) = intern_term(i,j) + exp (-((domain_len(t)^2*(i-cells(k,1))^2)+(j - cells(k,2))^2));
            end
      
            chemo_new(i,j) = dt * (D*((1/domain_len(t)^2)* (chemo(i+1)-2*chemo(i,j)+chemo(i-1,j))/(di^2)...
                + (chemo(i,j+1)- 2* chemo(i,j)+chemo(i,j-1))/dj  ) -  (chemo(i,j)*lam/ (2*pi*R^2)) * intern_term(i,j) +  ...
                kai*chemo(i,j)*(1-chemo(i,j)) + domain_len_der(t)/domain_len(t) *chemo(i,j) ) + chemo(i,j);
        end

    end

    %%%%%%% move cells
    order = randperm(size(cells,1)); % random permutation, we will pick cells in this order
    angle = zeros (1,size(cells,1)); % initialise angle vector

    %%%%% random walk, FALSE
%     for i = 1 : size(cells,1)
%         random_vector = rand(2,1)-[0.5,0.5]';
%         unit_random_vector = random_vector/norm(random_vector);
%         angle(order(i)) = atan(unit_random_vector(2)/unit_random_vector(1)); % find an angle
%         cells(i,1) = cells(i,1) + cos(angle(order(i))); % update x coordinate
%         cells(i,2) = cells(i,2) + sin(angle(order(i))); % update y coordinate
%     end

    %%%% self-propelled model
    for i = 1 : size(cells,1)
        x_cor = round(cells(i,1));
        y_cor = round(cells(i,2));
        x_grad = gradient ([chemo(x_cor),chemo(x_cor+1)]);
        y_grad = gradient ([chemo(y_cor),chemo(y_cor+1)]);
        angle(order(i)) = atan(y_grad/x_grad); % find an angle due to
        %gradient
        %%% add noise
        random_vector = rand(2,1)-[0.5,0.5]'; % generate a random vector between -1/2 and 1/2
        unit_random_vector = random_vector/norm(random_vector); % normalise it
        
        angle(order(i)) = angle(order(i)) + 0.5*atan(unit_random_vector(2)/unit_random_vector(1)); % add noise
     
        cells(i,1) = cells(i,1) + cos(angle(order(i))); % update x coordinate
        cells(i,2) = cells(i,2) + sin(angle(order(i))); % update y coordinate
    end
     
    
    chemo = chemo_new; %%%% update chemoattractant values
    t = t +dt; % update time step, maybe for loop instead, if this, then would need to have a while loop
    
    hold on
    scatter(cells(:,1),cells(:,2),'b','filled')


end

%%% Heat map of the chemoattractant
     HeatMap(chemo')
     
     
function L = domain_len(t)

    % Parameters 
    L0 = 300; % initial length \nu m
    a = 0.08; %h^-1\nu m^-1
    ts = -16; % h
    Linf = 870; % \nu m 

    L = L0 * ((Linf*exp(a*(t-ts)*Linf))/(Linf-1+exp(a*(t-ts))*Linf) + 1 - ...
        (Linf*exp(a*(-ts)*Linf))/(Linf-1+exp(a*(-ts)*Linf)));
    
    L =100; % fixed domain
end


function L = domain_len_der(t)

    % Parameters 
    L0 = 300; % initial length \nu m
    a = 0.08; %h^-1\nu m^-1
    ts = -16; % h
    Linf = 870; % \nu m 

    L = L0 * ((Linf^2*a*exp(a*(t-ts)*Linf))/(Linf-1+exp(a*(t-ts))*Linf) - ...
       (Linf^2*a*exp(2*a*(t-ts)*Linf))/(Linf-1+exp(a*(t-ts))*Linf));
   
   L =0; %%%% no change
end
