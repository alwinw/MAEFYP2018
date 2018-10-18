clear
close all

%Constant strength doublet method:

%Define the free stream speed:
u_infty = 1;

%Aerofoil:
Aerofoil = 0012;

%chord:
c = 1;

%Thickness:
t = c*Aerofoil/100;

%Half the mesh resolution - 1:
Npanel = 90;

%This is to avoid any issues when one panel is taken out:
if Npanel < 8
    
    Npanel = 8;
    
end

%In angular terms,
delta_beta = pi/(Npanel);
beta = 0:delta_beta:2*pi;
%2*pi implies that top and bottom are being simulated.

%Corresponding x-coordinates:
x = c/2*(1 - cos(beta));
x_t = x(1:Npanel+1);
x_b = x(Npanel+1:2*Npanel);

%x_t is for the top side and goes from the leading edge to the trailing
%edge. x_b is for the bottom and goes from the trailing edge to the leading
%edge. This is because of the setup of the matrices defined at the bottom
%of page 282 in the Katz and Plotkin book which incorporate the Kutta
%condition into the simulation.

%Calculate the z-coordinates (top half only initially):
z_t = 5*t*(0.2969*sqrt(x_t/c) + (-0.1260)*(x_t/c) + (-0.3516)*(x_t/c).^2 + 0.2843*(x_t/c).^3 + (-0.1036)*(x_t/c).^4);
z_b = -(z_t(end:-1:2));

%Panel corner coordinates:
x = horzcat(x_b,x_t);
z = horzcat(z_b,z_t);
z(1) = 0;
z(end) = 0;

%Show the aerofoil:
figure
plot(x,z);

%Calculate theta_i:
dz = diff(z);
dx = diff(x);
theta = atan2(dz,dx);

figure
hold on

for al = 1:1:6
    
    %Define the angle of attack:
    alpha = (al-1)*pi/180;
    
    sin_theta = -sin(theta);
    cos_theta = cos(theta);
    %Definition of atan2 means angles have the wrong sign, hence the fix
    %above.
    
    %An alternative method would be to use the formulation provided by page
    %265 of Katz and Plotkin. However, it is easy to get signs of angles 
    %wrong with this method.
    
    %Note - those referring to the code in the back will find sin(alpha) is
    %missing a sign. Because of this, the matrices used for the
    %transformations to and from the panel coordinate system are slightly 
    %different to those defined by equations 11.23 and 11.23a. Again, it's 
    %to do with the definition of atan2.
    
    %Collocation points (loocated halfway between panel corners):
    x_c = dx/2 + x(1:end-1);
    z_c = dz/2 + z(1:end-1);
    
    %Check the normal vector is correctly calculated:
    % hold on
    % plot(x_c+sin_theta,z_c+cos_theta,'x')
    % plot(x_c,z_c,'x')
    
    %By equation 11.3, the normal unit vector is defined as follows:
    n_hat = [sin_theta;cos_theta];
    %By equation 11.3a, the tangential unit vector is defined as follows:
    t_hat = [cos_theta;-sin_theta];
    
    %Number of collocation points:
    Ni = length(x_c);
    
    %Number of source points:
    Nj = length(x) - 1;
    
    %Create arrays with space for Kutta condition results:
    a = zeros(Ni+1,Nj+1);
    b = a;
    RHS = zeros(Ni+1,1);
    
    %Set the bottom row of the 'a' matrix to incorporate the Kutta
    %condition:
    a(Ni+1,1) = 1;
    a(Ni+1,Nj) = -1;
    a(Ni+1,Nj+1) = 1;
    
    %Complete the job for the RHS:
    RHS(Ni+1,1) = 0;
    
    %Rather than calculate X2 time and time again, it is in fact possible 
    %to calculate X2 simply by
    X2 = sqrt(dx.^2 + dz.^2);
    %This is because X2 is in the panel coordinate system.

    %Relocation of collocation points to inside the body:
    CP = 0;
    
    if CP == 1
        
        x_c = x_c - sin_theta.*X2*0.0125;
        z_c = z_c - cos_theta.*X2*0.0125;

    end
    
    for i = 1:Ni
        
        for j = 1:Nj
            
            %Equation 11.23a:
            Co = [cos_theta(j),-sin_theta(j);sin_theta(j),cos_theta(j)]*...
                [x_c(i) - x(j);z_c(i) - z(j)];
            
            %Radii:
            R1 = sqrt(sum(Co.^2));
            R2 = sqrt((Co(1) - X2(j))^2 + Co(2)^2);
            
            if i == j && CP == 0
                
                %Equations 11.34 and 11.35:
                up = 0;
                wp = -1/(pi*Co(1));
                
            else
                
                %Equations 11.30 and 11.31:
                up = 1/(2*pi)*(Co(2)/(R1^2) - Co(2)/(R2^2));
                wp = -1/(2*pi)*(Co(1)/(R1^2) - (Co(1)-X2(j))/(R2^2));
                
            end
            
            %Equation 11.23:
            Vec = [cos_theta(j),sin_theta(j);-sin_theta(j),cos_theta(j)]*...
                [up;wp];
            
            %By equation 11.27, the influence coefficient is as follows:
            a(i,j) = transpose(Vec)*n_hat(:,i);
            
            %Equation 11.37 (without the true values of mu included):
            b(i,j) = transpose(Vec)*t_hat(:,i);
            
            %According to equation 11.6b, the RHS is calculated as follows:
            RHS(i) = -u_infty*[cos(alpha),sin(alpha)]*n_hat(:,i);
            
        end
        
        %Influence of the wake panel:
        r_squared = (x_c(i)-x(end))^2 + (z_c(i) - z(end))^2;
        
        %For the wake panel, x_2 or x(Nj + 1) is equal to infinity. Also,
        %z(Nj + 1) is equal to zero. x_1 is simply x(Nj) and z_2 is zero.
        %Hence, by equations 11.30 and 11.31, the velocity components are
        Vec = 1/(2*pi*r_squared)*[z_c(i);-(x_c(i)-x(end))];
        %No conversion between coordinate systems is required as by
        %equations 11.23 and 11.23a, if alpha_i = 0 (which it will in this
        %case since dz = 0, then u = u_p and w = w_p.
        
        %Note - the second terms in both 11.30 and 11.31 go to zero by
        %virtue of the fact that the wake panel endpoint (x_2) is at
        %infinity.
        
        %Complete the 'a' matrix:
        a(i,Nj+1) = transpose(Vec)*n_hat(:,i);
        
        %Equation 11.37 (without the true values of mu included):
        b(i,Nj+1) = transpose(Vec)*t_hat(:,i);
        
    end
    
    %Remove one of the panels to help with singular matrices. This is a 
    %known problem when applying the Neumann boundary condition. Since all 
    %the other panels are trying to prevent flow in/out of the aerofoil, 
    %mass continuity should mean there is also negligible flow through the 
    %missing panel.    
    a_test = a;
    b_1 = horzcat(b(:,1:Npanel-4),b(:,Npanel-2:end));
    b = vertcat(b_1(1:Npanel-4,:),b_1(Npanel-2:end,:));
    a_1 = horzcat(a(:,1:Npanel-4),a(:,Npanel-2:end));
    a = vertcat(a_1(1:Npanel-4,:),a_1(Npanel-2:end,:));
    RHS = vertcat(RHS(1:Npanel-4),RHS(Npanel-2:end));
    n_hat = horzcat(n_hat(:,1:Npanel-4),n_hat(:,Npanel-2:end));
    t_hat = horzcat(t_hat(:,1:Npanel-4),t_hat(:,Npanel-2:end));
    x_c = horzcat(x_c(1:Npanel-4),x_c(Npanel-2:end));
    z_c = horzcat(z_c(1:Npanel-4),z_c(Npanel-2:end));
    cos_theta = horzcat(cos_theta(1:Npanel-4),cos_theta(Npanel-2:end));
    sin_theta = horzcat(sin_theta(1:Npanel-4),sin_theta(Npanel-2:end));
    
    %Solve the linear set of equations:
    mu = a\RHS;
    
    %Once the strengths of the source elements has been determined, the u
    %and w velocity components can be calculated:
    %Q_ti = zeros(Ni,1);
    Q_ti = zeros(Ni-1,1);
    
    %Equation 11.37:
    q_ti = b*mu;
    
    %Since the number of panels has been reduced by 1, the for loop only
    %goes to Ni - 1, with the final point being shifted to Ni-1 (and by 
    %extension the Ni - 1 term ----> Ni -2).
    
    %for i = 1:Ni
    for i = 1:Ni-1
        
        %Equation 11.37:
        temp = q_ti(i);
        
        %Add in the effect of panel on itself (equation 11.38):
        if i == 1
            
            %Calculate l:
            l = sqrt((x_c(2) - x_c(1))^2 + (z_c(2) - z_c(1))^2);
            
            %Calculate Vp:
            Vp = 0.5*(mu(2) - mu(1))/l;
            
        %elseif i == Ni
        elseif i == Ni-1
            
            %Calculate l:
            %l = sqrt((x_c(Ni) - x_c(Ni-1))^2 + (z_c(Ni) - z_c(Ni-1))^2);
            l = sqrt((x_c(Ni-1) - x_c(Ni-2))^2 + (z_c(Ni-1) - z_c(Ni-2))^2);
            
            %Calculate Vp:
            %Vp = 0.5*(mu(Ni) - mu(Ni-1))/l;
            Vp = 0.5*(mu(Ni-1) - mu(Ni-2))/l;            
            
        else
            
            %Calculate l:
            l = sqrt((x_c(i+1) - x_c(i-1))^2 + (z_c(i+1) - z_c(i-1))^2);
            
            %Calculate Vp:
            Vp = 0.5*(mu(i+1) - mu(i-1))/l;            
            
        end
        
        Q_ti(i) = cos(alpha)*cos_theta(i)*u_infty - ...
            sin(alpha)*sin_theta(i)*u_infty + temp + Vp;
        
    end
    
    %Pressure coefficient:
    Cp = 1 - Q_ti.^2/u_infty^2;
    
    %plot(x_c/c,Cp,'k',x,z,'--k');
    plot(x_c/c,Cp,x/c,z,'--k');
    
end

xlabel('x/c');
ylabel('Pressure coefficient c_P');
set(gca,'Ydir','reverse')
%ylim([-5 1.2])