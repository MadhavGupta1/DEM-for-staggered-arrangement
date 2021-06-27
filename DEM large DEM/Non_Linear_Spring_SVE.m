% MATLAB codes for Finite Element Analysis
% It consist of an Algorithm which consist Displacement Increment
clear
clear all
close all
format long g

%Structure p
p=struct();

%% Probem Defination

% Number of Tablets in x-direction
Nx = 3;
% Number of Tablets in y-direction
Ny = 4;

% Properties for Shear Stress-Strain Curve of each spring
G=100*10^6;
H=8*10^6;
tau_y=16*10^6;
gamma_u=3;
gamma_p_s=2; 
phi=0.9;
tau_s = tau_y + H*gamma_p_s;

% Parameters of the Tablets
rho_mean=4;
k_mean=0.5;
t=1;

ratio = 0.05;
delta_rho = ratio*rho_mean;

% Normal Distribution to vary center position of Tablet
z = normrnd(0,sqrt(2)*delta_rho,[1,(Nx+3)*Ny]);    %sqrt(2)* is removed

% Finding Even Number in Ny array 
j = 1:1:Ny;
iseven_2=rem(j,2);

% d: offset to the tablet given to each rho 
d = iseven_2*(1-k_mean)*rho_mean;

% d: For columnar Tablets (put inside the loop)
% d = iseven_2(j)*(rand*rho_mean);

% Total Number of Tablets
Number_of_Tablets = Nx*Ny + sum(iseven_2(:) == 0);
% Total Number of Non-Linear springsn
Number_of_Springs = Nx*(Ny-1)*2;
% Total Number of Tablets in which Spring is Present
iseven_1 = iseven_2;
iseven_1(end) = [];
Number_of_Spring_Tablet = Nx*(Ny-1) + sum(iseven_1(:) == 0);

% Positions of Tablet Nodes
x = zeros(Nx+3,Ny);
y = zeros(Nx+3,Ny);
q=1;
for i=1:1:Nx+3
    for j=1:1:Ny 
        x(i,j)= d(j) + (i-1)*rho_mean*t + z(1,q)*t;
        y(i,j)= -t*(j-1);
        q=q+1;
    end
end

% rho: Aspect Ratio of each Tablet
rho = zeros(Nx+1,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            rho(i-1,j) = (x(i+1,j)-x(i-1,j))/(2*t);
        end 
    else
        for i=2:1:Nx+1
            rho(i-1,j) = (x(i+1,j)-x(i-1,j))/(2*t); 
        end   
    end
end

% x_corner: Defining Left Bottom Corner of each Tablet
x_corner = zeros(Nx+2,Ny);
y_corner = zeros(Nx+2,Ny);
for j=1:1:Ny 
    if iseven_2(j)==0
        x_corner(2,j) = x(2,j)-rho(1,j)*t/2;
        for i=2:1:Nx+2
            if i==2
                x_corner(i,j) = x_corner(2,j);
                y_corner(i,j) = y(i,j)-t/2;
            else
                x_corner(i,j) = x_corner(i-1,j) + rho(i-2,j)*t;
                y_corner(i,j) = y(i,j)-t/2;
            end
        end 
    else
        x_corner(2,j) = x(2,j)-rho(1,j)*t/2;
        for i=2:1:Nx+1
            if i==2
                x_corner(i,j) = x_corner(2,j);
                y_corner(i,j) = y(i,j)-t/2;
            else             
                x_corner(i,j) = x_corner(i-1,j) + rho(i-2,j)*t;
                y_corner(i,j) = y(i,j)-t/2;
            end
        end   
    end
end

% k: overlap Ratio of Each Interface
k = zeros(Number_of_Spring_Tablet,1);
n=1;
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2            
            if i==2
                k(n,2) = (x_corner(i+1,j)-x_corner(i,j+1))/(rho(i-1,j)*t); 
                n=n+1;               
            elseif i==Nx+2
                k(n,1) = (x_corner(i-1,j+1)+rho(i-2,j+1)*t-x_corner(i,j))/(rho(i-1,j)*t); 
                n=n+1;                   
            else
                k(n,1) = (-x_corner(i,j)+x_corner(i,j+1))/(rho(i-1,j)*t); 
                k(n,2) = (x_corner(i+1,j)-x_corner(i,j+1))/(rho(i-1,j)*t);
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            if i==Nx+1
                k(n,1) = (x_corner(i+1,j+1)-x_corner(i,j))/(rho(i-1,j)*t);
                k(n,2) = (x_corner(i,j)+rho(i-1,j)*t-x_corner(i+1,j+1))/(rho(i-1,j)*t); 
                n=n+1;               
            else
                k(n,1) = (x_corner(i+1,j+1)-x_corner(i,j))/(rho(i-1,j)*t); 
                k(n,2) = (x_corner(i+1,j)-x_corner(i+1,j+1))/(rho(i-1,j)*t);
                n=n+1;
            end
        end        
    end
end

% Node_coordinate: Coordinates of each Tablet
Node_coordinate = zeros(Number_of_Tablets,2);
n=1;
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            Node_coordinate(n,:) = [x(i,j) y(i,j)];
            n=n+1;
        end 
    else
        for i=2:1:Nx+1
            Node_coordinate(n,:) = [x(i,j) y(i,j)];
            n=n+1;
        end   
    end
end

% Node_Number: number corresponding to each Tablet Node for Defining
%              Element Connections
Node_Number = zeros(Nx+2,Ny);
for j=1:1:Ny
    if iseven_2(j)==0
        for i=2:1:Nx+2
            if j==1
                Node_Number(i,j) = (i-1);
            else      
                Node_Number(i,j) = (i-1) + (Node_Number(Nx+1,j-1));
            end
        end 
    else
        for i=2:1:Nx+1
            if j==1
                Node_Number(i,j) = (i-1);
            else      
                Node_Number(i,j) = (i-1) + (Node_Number(Nx+2,j-1));
            end
        end   
    end
end

% Element_Nodes: Nodes cooresponding to each spring Element
Element_Nodes = zeros(Number_of_Springs,2);
% Element_Length: Interface Length for each Spring
Element_Length = zeros(Number_of_Springs,1);
n=1;
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2            
            if i==2
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];  
                Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i-1,j)*t;
                n=n+1;
            elseif i==Nx+2
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i-1,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i-1,j)*t;
                n=n+1;
            else
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i-1,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i-1,j)*t;
                n=n+1;
                Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];
                Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i-1,j)*t;
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i,j+1)];
            Element_Length(n,1) = k(Node_Number(i,j),1)*rho(i-1,j)*t;
            n=n+1;
            Element_Nodes(n,:) = [Node_Number(i,j) Node_Number(i+1,j+1)];
            Element_Length(n,1) = k(Node_Number(i,j),2)*rho(i-1,j)*t;
            n=n+1;
        end   
    end
end

%Prescribed_Dof: Nodes at which Displacement Boundary Condition is defined
n=1;
m=1;
for j=1:1:Ny
    if iseven_2(j)==0
        Prescribed_Dof(n,1) = Node_Number(2,j);
        n=n+1;
        Prescribed_Dof(n,1) = Node_Number(Nx+2,j);
        n=n+1;
        Disp_Dof(m,1) = Node_Number(Nx+2,j);
        m=m+1;
    end
end

%% Parameters for Stress-Strain Curve

% Properties of Stress-Strain curve
E = (phi^2/(1-phi))*(Element_Length.^2*G/t^2);
E_T_1 = (phi^2/(1-phi))*(Element_Length.^2*G/t^2)/(G/H + 1);
E_T_2 = (phi^2/(1-phi))*(Element_Length.^2*G/t^2)/((G/tau_s)*(gamma_p_s-gamma_u) + 1);
H_1 = E.*E_T_1./(E-E_T_1);
H_2 = E.*E_T_2./(E-E_T_2);
sigma_y_1 = tau_y*phi*Element_Length/t;
epsilon_p_u = (1./(Element_Length/t))*((1-phi)/phi)*(gamma_p_s); 

%% 

% Node coordinates
p.Node_coordinate=Node_coordinate;
% Num_Nodes: number of nodes 
p.Num_Nodes=size(p.Node_coordinate,1);

% Element_Nodes: connections at elements 
p.Element_Nodes=Element_Nodes;
% Num_Elements: number of Elements 
p.Num_Elements=size(p.Element_Nodes,1); 
% Element_Length: Length of each spring Element
p.Element_Length = Element_Length;
% Element_Area: Area of each spring Element
p.Element_Area = t*1*ones(p.Num_Elements);

% Element_Stiffness: truss stiffness values
p.Element_E = E;
% NDof: Degree of Freedom at each node
p.NDof = 1;
% GDof: total number of degrees of freedom
p.GDof=p.Num_Nodes;

%% Given Boundary Conditions

% Boundary condition for Homogeneous Solution
% Prescribed_Dof: Value of displacement is given 
p.Prescribed_Dof=Prescribed_Dof;
% Displacement Applied at Dof
p.Displacement_Dof = Disp_Dof; 
 
%% Initializing Matrix 

% Displacement: Displacement vector
p.Displacement=zeros(p.GDof,1);
% Displacment Increment:
Displacement_increment=zeros(p.GDof,1);
% Force: Force vector 
p.Force=zeros(p.GDof,1);
% Stiffness: Stiffness Matrix 
p.Stiffness=zeros(p.GDof,p.GDof);

% Number of Displacements Steps
n=1000000;

% u: Displacement for all Steps
u = zeros(p.Num_Nodes,n+1);
% Force: Force for all Steps
Force = zeros(p.Num_Nodes,n+1);
% Strain: for each Element for all Steps
Strain = zeros(p.Num_Elements,n+1);
% Stress: for each Element for all Steps
Stress = zeros(p.Num_Elements,n+1);

% Yeild Stress: at each Step
sigma_yeild = zeros(p.Num_Elements,1);
% Trial Stress: at each step
tr_sigma = zeros(p.Num_Elements,1);
% Plastic Strain: for all Steps
Plastic_Strain = zeros(p.Num_Elements,n+1);
% Change in Plastic Strain at each Step
Delta_plastic_strain = zeros(p.Num_Elements,1);

% Active Dof: Nodes at which Force is supposed to be zero
% As Displacement Increment is considered 
active_Dof = setdiff([1:p.GDof]', [p.Prescribed_Dof]);

%% Solving Spring Problem 

% Displacement to Strain Transformation Matrix
B = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    % Nodes of each Element 
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % Evaluation of Transformation matrix
    Le =(node_xx(1,2)-node_xx(1,1));
    
    % Length of each Element
    temp=[-sign(Le)/p.Element_Length(e) sign(Le)/p.Element_Length(e)];         
    B(e, el_node(1)) = B(e, el_node(1)) + temp(1);
    B(e, el_node(2)) = B(e, el_node(2)) + temp(2);
end

%Stress to Force Transformation at each Node
B_T = zeros(p.Num_Elements,p.Num_Nodes);       
for e=1:p.Num_Elements
    %Nodes of each Element
    el_node = p.Element_Nodes(e, 1:2);
    %x-coordinate of Node
    node_xx = p.Node_coordinate(el_node);
    % Evaluation of Transformation matrix
    Le = (node_xx(1,2)-node_xx(1,1));

    % Evaluation of Transformation matrix
    temp=[-sign(Le)*p.Element_Area(e)  sign(Le)*p.Element_Area(e)];         
    B_T(e, el_node(1)) = B_T(e, el_node(1)) + temp(1);
    B_T(e, el_node(2)) = B_T(e, el_node(2)) + temp(2);
end

% Loop Solving for Non-Linear Elastic Analysis

% Value for Accuracy
tol = 1.0e-6;
MAX_ITER = 1000;
trial_num = 0;
max_iter = MAX_ITER;
b_increment = 0.001;
b = b_increment;         %Increment
i=2;

% Cell array of colours 
C = {[0 0 0],[0.8 0.8 0.8],[1 1 0],[1 0 0]}; 
% Black:    'c'         : Elastic Region
% Grey:     'g'         : Elastic Unloading Region
% Yellow:   'y'         : Strain Hardening Region
% Red:      'r'         : Strain Softening Region
Colour(:,1) = C(1);
z=1;  

while i<n+1  

    % Displacement increment at each Step
    Displacement_increment(p.Displacement_Dof) = b;
    p.Displacement = Displacement_increment;

    % Stiffness Evaluation
    E_temp = p.Element_E;
    Stiffness = Global_Stiff_Matrix(p);
    % u_increment = p.Displacement;
    u_increment = Displacement(p,Stiffness);
    u(:,i) = u(:,i-1) + u_increment(:,1);
    % Force Evaluation due to Displacement Increment
    Force(:,i) = Force(:,i-1) + Stiffness*u_increment(:,1);
    
    %Displacement change after each step 
    Delta_u(:,1) = u(:,i)-u(:,i-1);
    %Displacement change after each iteration
    %delta_u(:,1) = zeros(p.GDof,1);
    
    % Convergence loop
    %conv = 10; 
    %iter = 0;
    
   % while conv > tol && iter < max_iter 

        % Evaluating new Global Stiffness Matrix
       % K_t = Global_Stiff_Matrix(p);
        
        % Displacement increment
       % Delta_u(:,1) = Delta_u(:,1) + delta_u(:,1);
        
        % Strain Increment
        % B=Disp_Strain(p,u(:,i));
        Delta_Strain(:,1) = B*Delta_u(:,1);
        
        I3 = find(Plastic_Strain(:,i-1)>=epsilon_p_u);
        I8 = setdiff([1:p.Num_Elements]', I3);
        
%         if ((isempty(I3)==0)&&(iter==0))
%             trial_num = trial_num + 1;
%         end
        
        % Trial Stress Evaluation
        tr_sigma(I8,1) = Stress(I8,i-1) + E(I8).*Delta_Strain(I8,1);
        tr_sigma(I3,1) = Stress(I3,i-1) + E(I3).*Delta_Strain(I3,1);
        
        % Yeild Stress Evaluation
        % For Hardening Region
        sigma_yeild(I8,i) = sigma_y_1(I8) + H_1(I8).*Plastic_Strain(I8,i-1);
       
        % For Softening Region
        sigma_yeild(I3,i) = sigma_y_1(I3) + H_1(I3).*epsilon_p_u(I3) + H_2(I3).*(Plastic_Strain(I3,i-1)-epsilon_p_u(I3));
  
        % Evaluation of Change in Plastic Strain 
        fr = abs(tr_sigma(:,1))-sigma_yeild(:,i);
        I0 = find(fr<=0);   
        I1 = setdiff([1:p.Num_Elements]', I0);
        
        Delta_plastic_strain(I8,1) = fr(I8)./(E(I8) + H_1(I8));
        Delta_plastic_strain(I3,1) = fr(I3)./(E(I3) + H_2(I3));
        Delta_plastic_strain(Delta_plastic_strain<0)=0;
        
        % Evaluation of Stress and Plastic Strain
        Plastic_Strain(:,i) = Plastic_Strain(:,i-1) + Delta_plastic_strain(:,1);
        Stress(I0,i) = tr_sigma(I0,1);
        Stress(I1,i) = tr_sigma(I1,1) - sign(tr_sigma(I1,1)).*E(I1).*Delta_plastic_strain(I1,1);
        
        % Updating Stiffness values According to different regions
        p.Element_E(I0) = E(I0);
        p.Element_E(I1) = E_T_1(I1);
        p.Element_E(I3) = E_T_2(I3);
        
        % Storing Colour information for Stress-Strain Curve
        I6 = I0;
        I4 = find(Plastic_Strain(:,i)>0);
        Common_Elements = intersect(I6,I4);
        I6 = setxor(I6,Common_Elements);
        Colour(I0,i) = C(2);
        Colour(I6,i) = C(1);
        Colour(I1,i) = C(3);
        Colour(I3,i) = C(4);   
        
        % Ressidual Evaluation
        % B_T = Stress_Force(p,u(:,i));
        %R = Force(:,i)-B_T'*Stress(:,i);

        % convergence is checked at active Dof 
        %conv = (norm(R(active_Dof)))^2;      
        
        % Displacement increment in each Iteration at Active Dof
        %delta_u(active_Dof,1) = K_t(active_Dof,active_Dof)\R(active_Dof);
        
        %iter = iter + 1; 
 %   end
%     
%     if trial_num==1
%         b=b/10; 
%         p.Element_E = E_temp;
%         trial_num = trial_num + 1;        
%         i=i-1; 
%     elseif (max_iter<(MAX_ITER+2*1000))&&(iter==max_iter||isnan(conv)||isinf(conv))
%         b=b/2;
%         max_iter = max_iter + 1000;
%         p.Element_E = E_temp;
%         i=i-1;        
%     else
%         max_iter=MAX_ITER;
        % Displacement Evaluation
        u(:,i) = u(:,i-1) + Delta_u(:,1);
        % Strain Evaluation
        % B=Disp_Strain(p,u(:,i));
        Strain(:,i) = B*u(:,i);
        
        if (any(Strain(:,i)<=0))&&(any(isnan(Stress(:,i)))==1||any(isinf(Stress(:,i)))==1)
            disp('Strain became negative')
            % Delete Elements Presents after Crack Formation
            Stress(:,i+1:n+1)=[];
            Strain(:,i+1:n+1)=[];
            u(:,i+1:n+1)=[];        %Displacement at which crack happens
            Force(:,i+1:n+1)=[];    %Force at the time of crack 
            break;
        end
    
        % The Crack has been taken place  if the below conditions are satisfied 
        if any(Stress(:,i) <= 0) %&& any(Strain(:,i) >= epsilon_p_max)
            disp('Crack is happenend in the Structure')
            % Delete Elements Presents after Crack Formation
            Stress(:,i+1:n+1)=[];
            Strain(:,i+1:n+1)=[];
            u(:,i+1:n+1)=[];        %Displacement at which crack happens
            Force(:,i+1:n+1)=[];    %Force at the time of crack
            break;
        end  
%     end
    i=i+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Formation of Stress-Strain Curve for whole Structure
%%%
p.Fixed_Dof = setdiff(p.Prescribed_Dof, p.Displacement_Dof);
Combined_L = zeros(size(p.Displacement_Dof,1),1);
q=1;
for j=1:1:Ny
    if iseven_2(j)==0
            Combined_L(q,1) = x(Nx+2,j)-x(2,j);
            q=q+1;
    end
end
Strain_combined = (u(p.Displacement_Dof,:)-u(p.Fixed_Dof,:))./(Combined_L);
if size(Strain_combined,1)==1
    Strain_combined_1 = Strain_combined;
else
    Strain_combined_1 = sum(Strain_combined)/size(Strain_combined,1);
end

total_t=0;
for j=1:1:Ny
    if j==1
        total_t = total_t + t/2;
    elseif j==Ny
        total_t = total_t + t/2;
    else
        total_t = total_t + t;
    end
end

Stress_combined = Force(p.Displacement_Dof,:)/(total_t*1);
if size(Stress_combined,1)==1
    Stress_combined_1 = Stress_combined;
else
    Stress_combined_1 = sum(Stress_combined);
end

figure(1)
plot(Strain_combined_1,Stress_combined_1,'k','LineWidth',1.5)
hold on
axis([0 inf 0 inf])
xlabel('Strain');
ylabel('Stress');
title('Stress-Strain Curve');
grid on;
grid minor;

% figure(2)
% plot(Strain(1,:),Stress(1,:),'r','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;
% 
% figure(3)
% plot(Strain(2,:),Stress(2,:),'b','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;
% 
% figure(4)
% plot(Strain(3,:),Stress(3,:),'b','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;
% 
% figure(5)
% plot(Strain(4,:),Stress(4,:),'b','LineWidth',1.5)
% hold on
% axis([0 inf 0 inf])
% xlabel('Strain');
% ylabel('Stress');
% title('Stress-Strain Curve for both the Elements');
% grid on;
% grid minor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tablet Formation of the Structure
%%%

figure(20)
for j=1:1:Ny 
    if iseven_2(j)==0
        for i=2:1:Nx+2
            hold on;
            if i==2
                L = rho(i-1,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            else
                L = rho(i-1,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor','b','EdgeColor','k','LineWidth',2);
            end
        end 
    else
        for i=2:1:Nx+1
            hold on;
            if i==2
                L = rho(i-1,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            else             
                L = rho(i-1,j)*t;
                rectangle('Position',[x_corner(i,j),y_corner(i,j),L,t],'FaceColor',[0 0 1],'EdgeColor','k','LineWidth',2);
            end
        end   
    end
end

n=1;
Col = size(Colour,2);
for j=1:1:Ny-1
    if iseven_2(j)==0
        for i=2:1:Nx+2 
            if i==2
                x=[x_corner(i,j)+rho(i-1,j)*t-Element_Length(n) , x_corner(i,j)+rho(i-1,j)*t];
                y=[y_corner(i,j) , y_corner(i,j)];
                line(x,y,'Color',cell2mat(Colour(n,Col)),'LineWidth',5)
                n=n+1;
            elseif i==Nx+2
                x=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y=[y_corner(i,j) , y_corner(i,j)];
                line(x,y,'Color',cell2mat(Colour(n,Col)),'LineWidth',5)
                n=n+1;
            else
                x=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
                y=[y_corner(i,j) , y_corner(i,j)];
                line(x,y,'Color',cell2mat(Colour(n,Col)),'LineWidth',5)
                n=n+1;
                x=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
                y=[y_corner(i,j) , y_corner(i,j)];
                line(x,y,'Color',cell2mat(Colour(n,Col)),'LineWidth',5)
                n=n+1;
            end        
        end
    else
        for i=2:1:Nx+1
            x=[x_corner(i,j) , x_corner(i,j)+Element_Length(n)];
            y=[y_corner(i,j) , y_corner(i,j)];
            line(x,y,'Color',cell2mat(Colour(n,Col)),'LineWidth',5)
            n=n+1;
            x=[x_corner(i,j)+Element_Length(n-1) , x_corner(i,j)+Element_Length(n-1)+Element_Length(n)];
            y=[y_corner(i,j) , y_corner(i,j)];
            line(x,y,'Color',cell2mat(Colour(n,Col)),'LineWidth',5)
            n=n+1;
        end   
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%