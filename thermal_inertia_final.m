%% Calculation of thermal inertia using MODIS images
% Xue and Varotsos (1995) and Chunfeng et al. (2013) proposed a thermal inertia model based on MODIS–LST data obtained four times a day and the time interval between measurements.This method has good potential for estimating soil thermal inertia in large areas. By simulating the daily changes in the surface temperature of MODIS and recovering the thermal inertia, it estimates the daily SHF
%The science of remote sensing has provided an in-situ method  to estimate the real thermal inertia of the soil. Xue and Varotsos  (1995), Chunfeng et al. (2013) and Ma (1997) derived a remote sensing model of the actual surface thermal inertia using a one- dimensional soil thermal conductivity equation and simplified  linear boundary conditions. A simple and operational thermal inertia model was developed using phase angle information and  the amount of daily temperature change. They used a second- order Fourier series approximation for LST to account for the effects of surface evaporation in regions with changes in surface moisture and vegetation cover.

%!! Important Note!! 
%In this code, 4 MODIS products are used for the daily temperature cycle model.
%which measure the land surface temperature at 10:30 AM/PM local solar time for the Terra satellite and 1:30 AM/PM for the Aqua satellite and are freely available as MOD11A1 and MYD11A1.
%To run this code, it is necessary to store the LST images of each of these
%four hours in a separate folder and number them in some way so that the image of each date can be read at the same time, for example, if we want to do modeling for a year. We have 365 images in each folder, which are numbered from 001 to 365.

%%!!Code execution steps
%1-Change the path of the LST and Albedo folder
%2-Enter the number of images you want in front of for.
%3-In this code, the LST images were multiplied by a scale factor of 0.02 to convert to temperature in Kelvin. If your images are in Kelvin, this part does not need to be executed.
%4-Enter the longitude and delta_UTC of the study area.
% 5-Enter the dimensions of the initial image in the reshape section
%6-Change the output path and save the image

% This code was written by Fahime Arabi. if you have any questions about it, I will answer you with the following email:
%Fahimearabi1993@gmail.com
%%
clc;
clear;
close all;
%Read Images from my computer
path_lst1='F:\Thesis\4-LST\6-DCY\YEARS\h1'  
path_lst10='F:\Thesis\4-LST\6-DCY\YEARS\h10'
path_lst13='F:\Thesis\4-LST\6-DCY\YEARS\h13'
path_lst22='F:\Thesis\4-LST\6-DCY\YEARS\h22'
path_albedo='F:\Thesis\5-G\modis\thermal inertia\albedo2020'

for i=1:365
    cd(path_lst1);
    a=dir('*.tif');
    namelst1=a(i).name;
    [img,Ref]=geotiffread([path_lst1, '\', namelst1]);
    strname=strsplit(namelst1,'_');
    ddatte=strname{1,3};
    info=geotiffinfo([path_lst1, '\', namelst1]);
    img=double(img);
    lst1=(img(:))';
    lst1=lst1*0.02;
    
    cd(path_lst10);
    a=dir('*.tif');
    namelst10=a(i).name;
    [img,Ref]=geotiffread([path_lst10, '\',  namelst10]);
    [B4,RefMatrx]=geotiffread( namelst10);
    InfoB4=geotiffinfo([ path_lst10, '\', namelst10]);
    img=double(img);
    lst10=(img(:))';
    lst10=lst10*0.02;
    
    cd(path_lst13);
    a=dir('*.tif');
    namelst13=a(i).name;
    [img,Ref]=geotiffread([path_lst13, '\', namelst13]);
    info=geotiffinfo([path_lst13, '\', namelst13]);
    img=double(img);
    lst13=(img(:))';
    lst13=lst13*0.02;
    
    cd(path_lst22);
    a=dir('*.tif');
    namelst22=a(i).name;
    [img,Ref]=geotiffread([path_lst22, '\', namelst22]);
    info=geotiffinfo([path_lst22, '\', namelst22]);
    img=double(img);
    lst22=(img(:))';
    lst22=lst22*0.02;
    
    cd(path_albedo);
    a=dir('*.tif');
    namealbedo=a(i).name;
    [img,Ref]=geotiffread([path_albedo, '\', namealbedo]);
    info=geotiffinfo([path_albedo, '\', namealbedo]);
    img=double(img);
    albedo=(img(:))';
   %%
    longitude=54.3569; %Enter the longitude of the study area
    delta_UTC=3.5; Enter %the delta_UTC of the study area
    d=i;
    LT=12;
    LSTM=15*delta_UTC;
    DC=2*pi/365;
    B=DC*(d+10)+0.033*sin(DC*(d-2));
    E0T=(9.87*(sin(2*B))+(7.6)*(sin(B-0.2)));
    TC=4*(longitude-LSTM)+E0T;
    LST=LT+(TC/60);
    tm=LST+1; %tm
    tm=13*3600;
    t2=13*3600;
    t1=1*3600;
    S0=1367; %S0 is solar constant (1367 W/m2)
    Ct=0.75; % CT is atmospheric transmission in the visible spectrum which is set to 

    b=(tan(0.0000727*(tm)))/(1-(tan(0.0000727*(tm))));
    omega1=atan((b)/(1+b));
    omega2=atan((b*(sqrt(2)))/(1+(b*(sqrt(2)))));
    t2=(13*3600);
    t1=(1*3600);
    d=i;
    Latitude=31.0833*0.0174532925 ;
    declination=23.45*sin((2*pi*(d-80))/(365));
    declination=declination*0.0174532925 ;
    phi=acos(tan(declination)*tan(Latitude));
    A1=((2/pi)*sin(declination)*sin(Latitude))+((1/2*pi)*cos(declination)*cos(Latitude)*(sin(2*phi)+(2*phi)));
    N=2
    A2=(((2*sin(declination)*sin(Latitude))/(N*pi))*(sin(N*phi)))+(((2*cos(declination)*cos(Latitude))/(pi*(N*N-1)))*(N*sin(N*phi)*cos(phi)-(cos(N*phi)*sin(phi))));
  p2=(A1*(cos((0.0000727*(t1))-omega1))-(cos((0.0000727*(t2))-omega1)))/sqrt((1+(1/b)+(1/(2*(b^2)))));
  A2=((sin(declination)*sin(Latitude)*sin(2*phi))/pi)+(((2*cos(declination)*cos(Latitude))/3*pi)*((2*sin(2*phi)*cos(phi))-(cos(2*phi)*sin(phi))));
  p3=(A2*((cos((0.0000727*(t2))-omega2))-(cos((0.0000727*(t1))-omega2))))/sqrt((2+((sqrt(2))/b)+(1/(2*(b^2)))))
  Coefficient=p2+p3;
 result = Coefficient;
   
 deltaT=(lst13-lst1); 
 a=albedo;
  pp=(1-a)*S0*0.75;
  pp2=deltaT*sqrt(0.0000727);
  ppp=pp./pp2;
  Thermal_inertia=ppp*Coefficient;

 %% Save the output
   Thermal_inertia=reshape(Thermal_inertia,160,181); %Enter the dimensions of the initial image 
    nnaammee=['Thermal_inertia_' num2str(ddatte) ];
   filenam = ['F:\Thesis\5-G\modis\thermal inertia\test\' nnaammee '.tif']  %Change the output path and save the image
   geotiffwrite( filenam,Thermal_inertia,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
   
end