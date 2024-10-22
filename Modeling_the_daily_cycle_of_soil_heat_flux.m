%% Modeling the daily cycle of soil heat flux using MODIS sensor images
% This code can estimate the soil heat flux image in 24 hours by having four soil heat flux images.
%Using SHF estimation by the energy balance equation and with the help of NDVI, LST, albedo, and emissivity products, SHF is estimated four times a day employing MODIS images.
% Then using these four SHF images, modeling is done to estimate the daily cycle.The code for estimating soil heat flux during the day and night has already been provided by the author.


% This code was written by Fahime Arabi. if you have any questions about it, I will answer you with the following email:
%Fahimearabi1993@gmail.com
%%
clc
clear all;
close all;
% Read Images from my computer
path_G1='F:\Thesis\5-G\modis\soil heat flux\2020\G-MYD\NIGHT'
path_G10='F:\Thesis\5-G\modis\soil heat flux\2020\G_MOD\1-DAY'
path_G13='F:\Thesis\5-G\modis\soil heat flux\2020\G-MYD\DAY'
path_G22='F:\Thesis\5-G\modis\soil heat flux\2020\G_MOD\2-NIGHT'


for i=1:365  %Enter the number i based on the number of input images.
    
    cd(path_G1);
    a=dir('*.tif');
    nameG1=a(i).name;
    [img,Ref]=geotiffread([path_G1, '\', nameG1]);
    strname=strsplit(nameG1,'_');
    ddatte=strname{1,2};
    info=geotiffinfo([path_G1, '\', nameG1]);
    img=double(img);
    G1=(img(:))';

    cd(path_G10);
    a=dir('*.tif');
    nameG10=a(i).name;
    [img,Ref]=geotiffread([path_G10, '\', nameG10]);
     [B4,RefMatrx]=geotiffread(nameG10);
    InfoB4=geotiffinfo([path_G10, '\', nameG10]);
    img=double(img);
    G10=(img(:))';
    
    cd(path_G13);
    a=dir('*.tif');
    nameG13=a(i).name;
    [img,Ref]=geotiffread([path_G13, '\', nameG13]);
    info=geotiffinfo([path_G13, '\', nameG13]);
    img=double(img);
    G13=(img(:))';
    
    cd(path_G22);
    a=dir('*.tif');
    nameG22=a(i).name;
    [img,Ref]=geotiffread([path_G22, '\', nameG22]);
    info=geotiffinfo([path_G22, '\', nameG22]);
    img=double(img);
    G22=(img(:))';
    
 data=cat(1,G10,G13,G22,G1);%Image Read
 %% In this section, sunrise, sunset and day length are calculated.
    d=i; %Enter the Latitude of your study area
    Latitude=31.0833;
    declination=23.45*sin((2*pi*(d-80))/(365));
 
    Axis=23.439*pi/180;
    j=pi/182.625;
    m=1-tan(Latitude*pi/180).*tan(Axis*cos(j*d));

     m(m>2)=2;
     m(m<0)=0;
     b=acos(1-m)/pi;
     hours=b*24;
     w=hours;  %Length of day
     sunrise=12-(w/2); %Sunrise
     sunset=sunrise+w;  %Sunset
     ts=sunset- sunrise;
%% In this section, local noon time is calculated
 longitude=54.3569; %Enter the longitude of your study area
    delta_UTC=3.5;
    d=i;
    LT=12;
    GM=15*delta_UTC;
    DC=2*pi/365;
    B=DC*(d+10)+0.033*sin(DC*(d-2));
    E0T=(9.87*(sin(2*B))+(7.6)*(sin(B-0.2)));
    TC=4*(longitude-GM)+E0T;
    G=LT+(TC/60);
    tm=G+1; %tm
    tm=tm- sunrise;
    u=(pi/w)*(ts-tm);
    k=(w/pi)*atan((pi/w)*(ts-tm))
%% calculation mean and amplitude 
 for i= 1:115116 %In the line below, enter the product of the number of rows and columns of your image against i.
    y=[data(:,i)];
    f=[1, cos(((pi/w)*((10- sunrise)-tm))); 1, cos(((pi/w)*((13- sunrise)-tm)));1, cos(((pi/w)*((22- sunrise)-tm)))*exp(-((22- sunrise)-ts)/k);1, cos(((pi/w)*((25- sunrise)-tm)))*exp(-((25- sunrise)-ts)/k)];;               
    a=(inv(transpose(f)*f))*transpose(f)*double(y);
    result(:,i)=a;
end
mean=result(1,:);
amplitude=result(2,:);
%%
Gmax=mean+amplitude*cos((pi/w)*(tm-tm));
Gsrd=mean+amplitude*cos((pi/w)*((sunrise-1)-(G+1)));
Gsrd1 =(mean)+((amplitude*cos((pi/w)*(sunset-(G+1))))*exp(-(24-sunset)/k));
Ga=(Gmax-Gsrd)/(cos(pi/4)+1);
G0=Gsrd+Ga.*cos(pi/4);
omrgaG=Ga.*(cos(u).*(Ga*cos(u)+G0-Gsrd1)+((pi/4).*sin(u).*(24-ts)*(G0-Gsrd1)))/(G0-Gsrd1-Ga.*((pi/4).*sin(u).*(24-ts)-cos(u)));
u=(pi/w)*(ts-tm);
u=(pi/w)*(ts-tm);
k=((Ga.*(cos(u))-omrgaG))/((Ga.*pi/w)*sin(u));
ts=round(ts)
tss=ts+1;
%% Estimation of daily soil heat flux cycle
for i=1:ts
    t=i;
    G =mean+amplitude.*cos((pi/w).*((t)-tm));
    G=reshape(G,318,362); % Enter the dimensions of your image
   nnaammee=['G_' num2str(ddatte) '_' num2str(i)];
   filenam = ['F:\Thesis\4-G\6-DCY\YEARS\hourly_size\' nnaammee '.tif']  % Change G folder paths.
   geotiffwrite( filenam,G,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
end
for i=tss:24
    t=i;
    G =(mean+omrgaG)+((amplitude*cos((pi/w)*(ts-tm))-omrgaG)*exp(-(t-ts)/k)); % Enter the dimensions of your image
    G=reshape(G,318,362); 
    nnaammee=['G_' num2str(ddatte) '_' num2str(i)];
   filenam = ['F:\Thesis\4-G\6-DCY\YEARS\hourly_size\' nnaammee '.tif']  %Change G folder paths.
   geotiffwrite( filenam,G,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
 
end
 end