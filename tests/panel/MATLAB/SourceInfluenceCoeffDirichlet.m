%% Influence coefficients:

function [b] = SIC(cos_theta,sin_theta,x_c,z_c,X2,x,z,CP)



Ni = length(x_c);
Nj = length(x) - 1;

%Create array:
b = zeros(Ni,Nj);

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

    %Radii squared:    
    R12 = Co_x.^2 + Co_z.^2; 
    R22 = (Co_x - X2).^2 + Co_z.^2;        
    
    %Angle parts of equation 11.64. Note eta angles are the theta angles of
    %figure 11.17.
    eta_1 = atan2(Co_z,Co_x);
    eta_2 = atan2(Co_z,Co_x - X2);
        
    %Equations 11.63:
    f = (Co_x.*log(R12) - (Co_x - X2).*log(R22) + 2*Co_z.*...
        (eta_2 - eta_1));
    b(i,:) = 1/(4*pi)*f;    
        
    %Throw in exception to handle cases where collocation point is on the
    %panel surface:
    if CP == 0
        
        b(i,i) = 1/(2*pi)*Co_x(i)*log(R12(i));
        
    end    
    
end