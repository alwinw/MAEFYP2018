function [x_c,z_c] = CollocationInternalise(X2,sin_theta,cos_theta,x_c,z_c,TestCp)

x1 = x_c;
z1 = z_c;

%Move collocation points inward by 0.05 panel lengths:
x_c = x_c - sin_theta*0.0125.*X2*4/4;
z_c = z_c - cos_theta*0.0125.*X2*4/4;

%Check the relocation of the collocation points:
plot(x_c,z_c,'x')

if TestCp == 1
    
    %Show that the angle is preserved:
    no = Npanel-5;
    testx = [x1(no),x_c(no),x_c(no)+sin_theta(no)];
    testz = [z1(no),z_c(no),z_c(no)+cos_theta(no)];

    figure
    plot(testx,testz,testx,testz,'x')

end