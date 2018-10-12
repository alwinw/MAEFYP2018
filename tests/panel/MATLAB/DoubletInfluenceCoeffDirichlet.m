%% Influence coefficients:

function [a] = DoubletInfluenceCoeffDirichlet(cos_theta,sin_theta,x_c,z_c,X2,x,z,CP)

Ni = length(x_c);
Nj = length(x) - 1;

%Create arrays:
c = zeros(Ni+1,Nj+1);
a = zeros(Ni,Nj);
%By equation 11.70, the Kutta condition can be included by modifying the
%doublet influence coefficients. This reduces the matrix size by one.

for i = 1:Ni

    %Establish collocation point in the coordinate systems of all the
    %panels using equation 11.23a of Katz and Plotkin:
    Co_x = cos_theta.*(x_c(i) - x(1:end-1)) - sin_theta.*(z_c(i) - z(1:end-1));
    Co_z = sin_theta.*(x_c(i) - x(1:end-1)) + cos_theta.*(z_c(i) - z(1:end-1));
    %Sin_theta and cos_theta are expressing the orientation of the panels
    %with respect to the body coordinate system. Note, in this work the
    %alpha_i of Katz and Plotkin (figure 11.17 and equation 11.23a) are
    %relabelled as theta_i in order to avoid confusion with angle of
    %attack. As a result, theta of figure 11.17 is also relabelled as eta.
    
    %Angle parts of equation 11.64. Note eta angles are the theta angles of
    %figure 11.17.
    eta_1 = atan2(Co_z,Co_x);
    eta_2 = atan2(Co_z,Co_x - X2);
                
    %Equation 11.64:    
    c(i,1:end-1) = -1/(2*pi)*(eta_2 - eta_1);
    
    %Throw in exception to handle cases where collocation point is on the
    %panel surface:
    if CP == 0
        
        c(i,i) = 0.5;
        
    end    
       
    %Influence of the wake panel (special case of 11.64):
    theta_1 = atan2((z_c(i)-z(end)),(x_c(i)-x(end)));
    
    %Modified from Katz and Plotkin. This is to be consistent with the use
    %of atan2.
    theta_2 = atan2(z_c(i)-z(end),x_c(i)-10000000*x(end));
    %Wake panel corner point not attached to trialing edge
    
    c(i,Nj+1) = 1/(2*pi)*(theta_1 - theta_2);
               
end

%Create a (using equation 11.70):
a(:,1) = c(1:end-1,1) - c(1:end-1,Nj+1);
a(:,end) = c(1:end-1,Nj) + c(1:end-1,Nj+1);
a(:,2:end-1) = c(1:end-1,2:end-2);
%This should yield a symmetric matrix which is easily invertible.