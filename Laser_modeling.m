%function laser_incomplete_MATLAB()

    %% %% Laser Heating of Block of Material %% %%
    clear all;
    close all;
    
    %% Assumptions
    % - Stress power is negligible i.e. sigma:grad(v) -> 0
    % - Deformations are negligible in the material derivative i.e.
    %   D(theta)/Dt -> thetadot
    % - Energy is only stored as heat (no mechanical work) i.e. all changes
    %   in internal energy result in temperature change
    % - Fourier's Law holds i.e. all heat flux comes from 
    %   internal conduction
    % - Constant density i.e. rho(t) = rho0
    % - Thermal conductivity is constant and isotropic
    
    %% Givens  
  
    % Intensive properties  
    rho0 = 6500;             % density in reference configuration, [kg/m3]
    C = 450;                % specific heat capacity, [J/K/kg]
    K = 135;                % isotropic thermal conductivity, [W/m/K]
    a = 25;                % absorption constant, [1/m]   
    
    % Extensive properties
    L = 0.05;                % cube dimension, [m]
    l = 0.01;                % laser dimension, [m]
    
    % Applied settings
    I0 = 2e8;               % initial laser intensity, [W/m2]
    thetaw = 300;           % fixed wall temperature, [K]
  
    % Initial Conditions  
    theta0 = 300;           % initial temperature, [K]
    
    % Time discretization
    T  = 2;               % final time for simulation, [s]
    dt = 1e-3;               % time step, [s]
    N  = T/dt;               % # of time points, [unitless]
    t  = linspace(0,T,N);               % time points, [s], 1xnt
    
    % Space discretization
    M = 21;                % # of nodes per dimension (9261 nodes total)
    h = L/(M-1);                % discritized element size, [m]
    disp(['h = ', num2str(h)]);
    x = 0:h:L;
    y = 0:h:L;
    z = 0:h:L;
   
    %% Construct Mask Arrays
    
    % Construct x-edge/y-edge mask for Dirichlet BCs                                                                           
    dirichlet_mask = zeros(M, M, M);
    dirichlet_mask(1, :, :) = 1; % left  edge (x-edge)                                                                   
    dirichlet_mask(M, :, :) = 1; % right edge (x-edge)                                                                   
    dirichlet_mask(:, 1, :) = 1; % back  edge (y-edge)                                                                   
    dirichlet_mask(:, M, :) = 1; % front edge (y-edge)   
 
    % Construct z-edge mask for Neumann BCs                                                                                    
    neumann_mask = zeros(M, M, M);
    neumann_mask(:, :, 1) = 1; % bottom edge (z-edge)                                                                          
    neumann_mask(:, :, M) = 1; % top edge (z-edge)                                                                             

    % Construct mask for obtaining nodes that neighbor Neumann nodes                                                           
    neumann_neighbors = zeros(M,M,M);
    neumann_neighbors(:, :, 2) = 1; % neighbor to a bottom edge                                                                
    neumann_neighbors(:, :, M-1) = 1; % neighbor to a top edge                                                                  

    % Define laser area                                                                                                        
    laser_nodes = int16(floor(l/h)); % number of nodes per dimension for the laser                                         
    laser_start = int16(floor((M - laser_nodes)/2));
    laser_end = laser_start + laser_nodes;
    fprintf('The value of laser nodes is: %.2f\n',laser_nodes );

    laser_mask = zeros(M,M,M);
    laser_mask(laser_start:laser_end, laser_start:laser_end, :) = 1;
  
    % Reshape masks for application to theta as a vector                                                                       
    dirichlet_mask    = reshape(dirichlet_mask, [M^3,1]);
    neumann_mask      = reshape(neumann_mask, [M^3,1]);
    neumann_neighbors = reshape(neumann_neighbors, [M^3,1]);
    laser_mask        = reshape(laser_mask, [M^3,1]);
    
    % Cast to logical arrays
    dirichlet_mask    = logical(dirichlet_mask);
    neumann_mask      = logical(neumann_mask);
    neumann_neighbors = logical(neumann_neighbors);
    laser_mask        = logical(laser_mask);
    
    %% Initialize storage variables
    Tavg_bot = ones(1, M)*theta0;  % avg. temperature of the bottom plane
    Tavg_top = ones(1, M)*theta0;  % avg. temperature of the top plane

    % Initialize the temperature in the block
    theta = theta0*ones(M^3,N); % M^3 x N
    
    % Define height at each node
    height = zeros(M^3,1);
    for n=1:length(height)
        height(n)=L-(n/(M^2+1))*h;
    end
    
    

    %% Construct Second Derivative Helper Matrix, D2                                                                         
    left_neighbors   = ones(M^3-1, 1);
    right_neighbors  = ones(M^3-1, 1);
    front_neighbors  = ones(M^3-M, 1);
    back_neighbors   = ones(M^3-M, 1);
    for i=1:M
        for j=1:M
            for k=1:M
  
              % Within this loop n represents a node number                                                                 
              n = i + (j-1)*M + (k-1)*M^2;
  
              if i == 1 % left edge (x-edge)                                                                               
                  if n-1 > 0
                      left_neighbors(n-1) = 0;
                  end
              end
  
              if j == 1 % back edge (y-edge)                                                                               
                  if n-M > 0
                      back_neighbors(n-M) = 0;
                  end
              end
            end
        end
    end
                      
    right_neighbors = left_neighbors;
    front_neighbors = back_neighbors;
  
    % Correct the D2 matrix to account for Dirichlet BCs                                                                      
    D2 = zeros(M^3, M^3);
    D2 = D2 + diag(-6*ones(1, M^3),   0);
    D2 = D2 + diag(left_neighbors,   -1);
    D2 = D2 + diag(right_neighbors,  1);
    D2 = D2 + diag(back_neighbors,   -M);
    D2 = D2 + diag(front_neighbors,  M);
    D2 = D2 + diag(ones(1, M^3-M^2), -M^2); % downstairs neighbors
    D2 = D2 + diag(ones(1, M^3-M^2), M^2); % upstairs neighbors
    D2 = sparse(D2); % convert to datatype sparse for speed
    
    % Uncomment the below only when M=4

    %{
    figure(1)
    spy(D2)
    xlabel('Nodes', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('Nodes', 'Interpreter', 'latex', 'FontSize', 20)
    title('Laplacian Matrix (4 x 4 x 4)', 'Interpreter', 'latex', 'FontSize', 20)
    saveas(gcf, 'spy.png')
    %}
    
    
    %% Time-stepping through the governing equation

    for n=1:length(t)-1
        
        % Within this loop n represents a time point
        disp(strcat(num2str(n), " of ", num2str(N)));

        % Get the values for the nodes with Neumann B.C.s
        neighbors = theta(neumann_neighbors, n);

        % Internal Nodes
        diffusion = (K/h^2)*D2*theta(:,n);
        beer_lambert = a*I0*exp(-a*height);
        theta(:,n+1) = theta(:,n) + dt/(rho0*C) * diffusion;
        theta(laser_mask,n+1) = theta(laser_mask,n+1) + dt/(rho0*C) * beer_lambert(laser_mask);
        

        % Dirichlet Fixed Temp. B.C.s at walls
        theta(dirichlet_mask, n+1) = thetaw;
        
        
        % Neumann Insulating B.C.s at top and bottom of block
        theta(neumann_mask, n+1) = neighbors;
                   
        % Store average temperatures of top and bottom planes
        Tavg_bot(n+1) = mean(theta(1:M^2, n+1),'all'); % z = 0
        Tavg_top(n+1) = mean(theta(M^3-(M-1)^2:M^3, n+1),'all'); % z = L 
        
    end
    
    % Final temp. profile of the top plane
    Tbot = reshape(theta(1:M^2,end), [M, M]); 
    % Final temp. profile of the top plane
    top_start = M^3 - (M^2-1);
    Ttop = reshape(theta(top_start:end,end),[M, M]); 
    
    %% DELIVERABLES %%
    
    %% Beer Lambert Intensity
    figure(2)
    plot(z, I0*exp(-a*z), 'LineWidth', 3)  
    title('Beer-Lambert Intensity','Interpreter', 'latex', 'FontSize', 20)
    xlabel('Depth (m)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('Laser Intensity ($W/m^2$)','Interpreter', 'latex', 'FontSize', 20)
    legend('I(z)','location','best', 'Interpreter', 'latex', 'FontSize', 14)
    saveas(gcf, 'beer_lambert.png')

    %% Beer Lambert volumetric power absorption
    figure(3) 
    plot(z, a*I0*exp(-a*z), 'LineWidth', 3)
    title('Beer-Lambert Volumetric Power Absorption','Interpreter', 'latex', 'FontSize', 20)
    xlabel('Depth (m)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('Volumetric Power Absorption ($W/m^3$)', 'Interpreter', 'latex', 'FontSize', 20)
    legend('$I_{abs}(z)$','location','best', 'Interpreter', 'latex', 'FontSize', 14)
    saveas(gcf, 'beer_lambert_vol.png')

    %% Average of top and bottom plane
    figure(4) 
    plot(t, Tavg_bot, t, Tavg_top, 'LineWidth', 3) 
    title('Avg. Temperatures of Top and Bottom Planes','Interpreter', 'latex', 'FontSize', 20)
    xlabel('Time (s)','Interpreter', 'latex', 'FontSize', 20)
    ylabel('Average Temperature (K)','Interpreter', 'latex', 'FontSize', 20)
    legend('bottom','top','location','best', 'Interpreter', 'latex', 'FontSize', 14)
    saveas(gcf, 'avg_temp.png')

    %% 2D grid of the temperatures at the top plane (final time)
    figure(5)
    [xq, yq] = meshgrid(x);
    surf(yq, xq, Ttop)
    title('Final Temperature Profile of Block Top Plane', 'Interpreter', 'latex', 'FontSize', 20)
    xlabel('x (m)', 'Interpreter', 'latex', 'FontSize', 20)
    ylabel('y (m)', 'Interpreter', 'latex', 'FontSize', 20)
    view(2)
    caxis([300, 1e3])  % the colorbar range from 300 to 1,000 K
    c = colorbar;
    colormap jet
    c.Label.String = 'Temperature (K)';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 14;
    saveas(gcf, 'top_plane_end_temp.png')

    %% 2D grid of the temperatures in the bottom plane (final time)
    figure(6)  
    surf(yq, xq, Tbot)
    title('Final Temperature Profile of Block Bottom Plane', 'Interpreter', 'latex', 'FontSize', 20)
    xlabel('x (m)','Interpreter', 'latex', 'FontSize', 20)
    ylabel('y (m)','Interpreter', 'latex', 'FontSize', 20)
    view(2)
    caxis([300,1e3])  % the colorbar range from 300 to 1,000 K
    c = colorbar;
    colormap jet
    c.Label.String = 'Temperature (K)';
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 14;
    saveas(gcf, 'bottom_plane_end_temp.png')


    %% Movie, for your amusement
%     xslice = x(ceil(M/2)); % location of y-z planes
 %    yslice = x(ceil(M/2)); % location of x-z plane
 %    zslice = [0, x(ceil(M/3)), x(ceil(2*M/3)), L]; % location of x-y planes
    
  %   v = VideoWriter('temperature_evolution', 'MPEG-4');
  %   open(v) 
   %  for n = 1:50:N
         
  %       disp(strcat(num2str(n), " of ", num2str(N)))
  %       figure(7)
  %       slice(x, y, z, reshape(theta(:,n), M, M, M),...
  %             xslice, yslice, zslice, 'cubic')  
  %      caxis([300,1e3])  % the colorbar range from 300 to 1,000 K
  %       c = colorbar;
    %     colormap jet
   %      c.Label.String = 'Temperature (K)';
  %%       c.Label.Interpreter = 'latex';
  %       c.Label.FontSize = 14;
    %     title('Temperature Evolution of Block Undergoing Laser Heating',...
   %            'Interpreter', 'latex', 'FontSize', 18)
   %      xlabel('x (m)', 'Interpreter', 'latex', 'FontSize', 20)
    %     ylabel('y (m)', 'Interpreter', 'latex', 'FontSize', 20)
    %     zlabel('z (m)', 'Interpreter', 'latex', 'FontSize', 20)
    %     hold off
   %      frame = getframe(figure(7));
   %      writeVideo(v,frame);
         
   %  end
  %   close(v);
