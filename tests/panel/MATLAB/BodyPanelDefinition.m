%% Body panel definitions - symmetric NACA aerofoils only

function [x,z,x_c,z_c,n_hat,t_hat,sin_theta,cos_theta,X2] = ...
    BodyPanelDefinition(c,t,Npanel,m,p)

%In angular terms,
delta_beta = pi/(Npanel);
beta = 0:delta_beta:2*pi;
%2*pi implies that top and bottom are being simulated.

%Corresponding x-coordinates:
x = c/2*(1 - cos(beta));

%Top part:
x_t = x(1:Npanel+1);

if m == 0 && p == 0

    x_b = x(Npanel+1:2*Npanel);
    
    %x_t is for the top side and goes from the leading edge to the trailing
    %edge. x_b is for the bottom and goes from the trailing edge to the 
    %leading edge. This is because of the setup of the matrices defined at 
    %the bottom of page 282 in the Katz and Plotkin book which incorporate 
    %the Kutta condition into the simulation.
    
    %Calculate the z-coordinates (top half only initially):
    z_t = 5*t*(0.2969*sqrt(x_t/c) + (-0.1260)*(x_t/c) + ...
        (-0.3516)*(x_t/c).^2 + 0.2843*(x_t/c).^3 + ...
        (-0.1036)*(x_t/c).^4);
    
    z_b = -(z_t(end:-1:2));
    
    %Panel corner coordinates:
    x = horzcat(x_b,x_t);
    z = horzcat(z_b,z_t);
    z(1) = 0;
    z(end) = 0;
    
else
    
    %Divide into two parts: that where x is less than p*c and that where is
    %equal to or greater than p*c. First step is determining transition:
    [~,x_transition] = min(abs(x_t - p*c));

    %This part is the same...
    z_t = 5*t*(0.2969*sqrt(x_t/c) + (-0.1260)*(x_t/c) + ...
        (-0.3516)*(x_t/c).^2 + 0.2843*(x_t/c).^3 + ...
        (-0.1036)*(x_t/c).^4);
    
    z_camb1 = m*x_t(1:x_transition-1)/p^2.*(2*p - ...
        x_t(1:x_transition-1)/c);
    grad = 2*m/(p^2)*(p - x_t(1:x_transition-1)/c);
    kappa1 = atan(grad);      
    
    z_camb2 = m*(c - x_t(x_transition:end))/((1 - p)^2).*(1 + ...
        x_t(x_transition:end)/c - 2*p);
    grad = 2*m/((1 - p)^2)*(p - x_t(x_transition:end)/c);
    kappa2 = atan(grad);
    
    kappa = horzcat(kappa1,kappa2);
    z_camb = horzcat(z_camb1,z_camb2);
    
    x_u = x_t - z_t.*sin(kappa);
    x_l = fliplr(x_t + z_t.*sin(kappa));
    
    z_u = z_camb + z_t.*cos(kappa);
    z_l = fliplr(z_camb - z_t.*cos(kappa));
    
    z = horzcat(z_l(1:end-1),z_u);
    x = horzcat(x_l(1:end-1),x_u);    
    
end

%Show the aerofoil:
figure
plot(x,z);

%Calculate theta_i:
dz = diff(z);
dx = diff(x);
theta = atan2(dz,dx);

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
hold on
plot(x_c+sin_theta,z_c+cos_theta,'x')
plot(x_c,z_c,'x')

%By equation 11.3, the normal unit vector is defined as follows:
n_hat = [sin_theta;cos_theta];
%By equation 11.3a, the tangential unit vector is defined as follows:
t_hat = [cos_theta;-sin_theta];

%Rather than calculate X2 time and time again, it is in fact possible to
%calculate X2 simply by
X2 = sqrt(dx.^2 + dz.^2);
%This is because X2 is in the panel coordinate system.