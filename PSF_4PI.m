% % % % % % % % %program body % % % % % % % % % % % % %
%Citations:
% https://opg.optica.org/josaa/fulltext.cfm?uri=josaa-9-12-2159&id=64091
% https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1959.0200
% https://opg.optica.org/josaa/fulltext.cfm?uri=josaa-27-10-2188&id=205731

% basic parameters
NA= 1.2; % numerical aperture of objective
n=1.33; % refractive index of immersion medium
% Individual PSF for multiplication
lambda785 = 790 * 10^-9;
lambda1029 = 1027 * 10^-9;
alpha=asin(NA/n); % maximum open angle of OL
k785=2*pi*n/lambda785; % wavenumber
k1029=2*pi*n/lambda1029;
li=sqrt(-1);

% image plane in Cartesian coordinate
% For displaying as 4Pi PSF figure. Keep as 2*10^-6 Lfocal
% For generating simulation PSF with 30nm steps, use 6*10^-6 Lfocal
Lfocal =1.5*10^-6; % observation dimension
Nx=100; % discretization of image plane
Ny=100; % discretization of image plane
x2=linspace(-Lfocal,Lfocal,Nx);
z2=linspace(-Lfocal,Lfocal,Ny);
[X2,Z2]= meshgrid(x2,z2);
Y2=0;

% polarization case
% ’1’: x-linear,’2’: y-linear,’3’: left circular
% ’4’: right circular,’5’: elliptical,’6’: radial,’7’: azimuthal
polar =1;
polarSTR=num2str(polar);

%normalization and steps of intergral
Ex2=0; % Ex-component in focal
Ey2=0; % Ey-component in focal
Ez2=0; % Ez-component in focal
Ex_2=0; % counter Ex-component in focal
Ey_2=0; % counter Ey-component in focal
Ez_2=0; % counter Ez-component in focal
Ex785det=0;
Ey785det=0;
Ez785det=0;
Ex785det2=0;
Ey785det2=0;
Ez785det2=0;

Ntheta =50;
Nphi =50;
deltatheta=alpha/Ntheta;
deltaphi=2*pi/Nphi;

% starting loop
for theta= 0*alpha:deltatheta:alpha %change starting theta to create annulus for Bessel beam, Gaussian beam should be a plane
    for phi= 0:deltaphi:2*pi
        
        % convertion funtion of polarization from ojective plane to imaging plane
        a=1+(cos(theta)-1)*(cos(phi))^2;
        b=(cos(theta)-1)*cos(phi)*sin(phi);
        c=sin(theta)*cos(phi);
        d=1+(cos(theta)-1)*(sin(phi))^2;
        e=-sin(theta)*sin(phi);
        ff=cos(theta);
        V=[a b c; b d e;-c -e ff];
        
        % incident beam polarization cases
        px =[1,0,1/sqrt(2),li/sqrt(2),2/sqrt(5),cos(phi),-sin(phi)];
        py =[0,1,li/sqrt(2),1/sqrt(2),li/sqrt(5),sin(phi),cos(phi)];
        pz=0;
        % selected incident beam polarziation
        P=[px(1,polar); py(1,polar); pz];
        % polarization in focal region
        PP=V*P;
        % numerical calculation of field distribution in focal region
        % Used to generate Illumination PSF        
        
        
        %Positive-Z propagation wave
        Ex2=Ex2+1i*sin(theta)*sqrt(cos(theta)).*PP(1,1).*exp(li*k1029*(Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ey2=Ey2+1i*sin(theta)*sqrt(cos(theta)).*PP(2,1).*exp(1i*k1029*(Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ez2=Ez2+1i*sin(theta)*sqrt(cos(theta)).*PP(3,1).*exp(1i*k1029*(Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        %Negative-Z or conuter-propagation wave
        Ex_2=Ex_2+1i*sin(theta)*sqrt(cos(theta)).*PP(1,1).*exp(li*k1029*(-Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ey_2=Ey_2+1i*sin(theta)*sqrt(cos(theta)).*PP(2,1).*exp(1i*k1029*(-Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ez_2=Ez_2+1i*sin(theta)*sqrt(cos(theta)).*PP(3,1).*exp(1i*k1029*(-Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;


        %Positive-Z propagation wave E1_785
        Ex785det=Ex785det+1i*sin(theta)*sqrt(cos(theta)).*PP(1,1).*exp(li*k785*(Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ey785det=Ey785det+1i*sin(theta)*sqrt(cos(theta)).*PP(2,1).*exp(1i*k785*(Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ez785det=Ez785det+1i*sin(theta)*sqrt(cos(theta)).*PP(3,1).*exp(1i*k785*(Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        
        % Negative-Z or counter-propagation wave E2_785
        Ex785det2=Ex785det2+1i*sin(theta)*sqrt(cos(theta)).*PP(1,1).*exp(li*k785*(-Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ey785det2=Ey785det2+1i*sin(theta)*sqrt(cos(theta)).*PP(2,1).*exp(1i*k785*(-Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
        Ez785det2=Ez785det2+1i*sin(theta)*sqrt(cos(theta)).*PP(3,1).*exp(1i*k785*(-Z2*cos(theta)+sin(theta).*(X2*cos(phi)+Y2*sin(phi))))*deltatheta*deltaphi;
    end
end


% Following how written into Properties of a 4Pi Confocal Fluorescence
% Microscope
% https://opg.optica.org/josaa/fulltext.cfm?uri=josaa-9-12-2159&id=64091
% Type A system with illumination PSF 785 and 1029 overlapped

% 4Pi-SRS system with multiplication of the interferring PSFs
|E1_1029 + E2_1029|.^2 .* |E1_785 + E2_785|.^2
Ix2=(abs(conj(Ex2+Ex_2).*(Ex2+Ex_2))) .* abs(conj(Ex785det+Ex785det2).*(Ex785det+Ex785det2));
Iy2=(abs(conj(Ey2+Ey_2).*(Ey2+Ey_2))) .* abs(conj(Ey785det+Ey785det2).*(Ey785det+Ey785det2));
Iz2=(abs(conj(Ez2+Ez_2).*(Ez2+Ez_2))) .* abs(conj(Ez785det+Ez785det2).*(Ez785det+Ez785det2));
I1=Ix2+Iy2+Iz2;

% Conventional SRS
% Ix2=abs(Ex2) .* abs(Ex785det);
% Iy2=abs(Ey2) .* abs(Ey785det);
% Iz2=abs(Ez2) .* abs(Ez785det);
% I1=Ix2+Iy2+Iz2;

% find maximum
MM1=max(max(I1));

%% plot xz view
hFig = figure(1);
set(hFig,'position',[300 300 300 300])
p=axes;
set(p,'Position',[0 0 1 1])
pcolor(X2*10^6,Z2*10^6,Ix2);
shading interp
title('PSF as a function of different polarization')
set(gca,'XTick',[-Lfocal 0 Lfocal]);
set(gca,'YTick',[-Lfocal 0 Lfocal]);
set(gca,'FontSize',12)
axis equal
axis tight
view(0,90);
hold on
colormap(hot)

%% Plotting axial profile
fig = figure('Units','pixels','Position',[1 1 512 550]);
I1 = I1 / max(I1(:));
hold on;
plot(x2*10^6,abs(I1(:,51)),'Color','red', 'LineWidth',1.5)
xlabel('z position / μm', 'FontName', 'Arial', 'FontSize', 28)
ylabel('Intensity / A.U.', 'FontName', 'Arial', 'FontSize', 28)
% % % % % % % % %end program% % % % % % % % % % % % % %
xlim([-1.3 1.3]);
ax = gca;
ax.FontSize = 25;
ax.FontName = 'Arial';
set(gca, 'LineWidth',1.1);
hold off;

