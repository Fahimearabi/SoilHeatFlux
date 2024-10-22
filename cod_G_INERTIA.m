%% Estimation of soil heat flux using thermal inertia
%Murray and Ver305 hoef (2007) proposed the harmonic analysis method, a physical model based on the harmonic analysis of soil surface temperature and independent of Rn. Their method only required the soil temperature of one layer instead of the temperature of multiple soil layers. 
% Therefore, the harmonic analysis method can be used not only for the point scale, but also for the regional scale. Determining errors and uncertainties in SHF estimation are more transparent and easier to interpret in the harmonic analysis method, because it is a physical method with semi empirical determination of inputs. 
% Note that this method can  only be used for the day. This analytical equation for estimating SHF, called the harmonic analysis method, is based on the harmonic analysis of the soil surface temperature.
% To run this code, images of land surface temperature at 1, 10, 13 and 22 hours are needed. Also needed is the image of the thermal inertia time series for which the code provided by the author is available.
% This code was written by Fahime Arabi. if you have any questions about it, I will answer you with the following email:
%Fahimearabi1993@gmail.com
%%
clc
close all
clear all
%Read Images from my computer
path_lst1='F:\Thesis\5-G\modis\lst\2020\2-LST_SSA_MYD\4-NIGHT_TIFF'
path_lst10='F:\Thesis\5-G\modis\lst\2020\1_LST_SSA_MOD\1-Day'
path_lst13='F:\Thesis\5-G\modis\lst\2020\2-LST_SSA_MYD\3-DAY_TIFF'
path_lst22='F:\Thesis\5-G\modis\lst\2020\1_LST_SSA_MOD\2-Night'
path_albedo='F:\Thesis\5-G\modis\albedo\image\albedo2020'
path_inertia='F:\Thesis\5-G\modis\DTG\G_INERTIA\Thermal_inertia'
%% 
for i=1:360 %Enter the number of images you want to process
    
    cd(path_lst1);
    a=dir('*.tif');
    namelst1=a(i).name;
    [img,Ref]=geotiffread([path_lst1, '\', namelst1]);
    strname=strsplit(namelst1,'_');
    ddatte=strname{1,2};
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
    lst10=lst10;
    
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
    lst22=lst22;
    
   data=cat(1,lst10,lst13,lst22,lst1);%Image Read
 
 cd(path_inertia);
    a=dir('*.tif');
    nameinertia=a(i).name;
    [img,Ref]=geotiffread([path_inertia, '\', nameinertia]);
    info=geotiffinfo([path_inertia, '\', nameinertia]);
    img=double(img);
    inertia=(img(:))';
    inertia=inertia;

    cd(path_albedo);
    a=dir('*.tif');
    namealbedo=a(i).name;
    [img,Ref]=geotiffread([path_albedo, '\', namealbedo]);
    info=geotiffinfo([path_albedo, '\', namealbedo]);
    img=double(img);
    albedo=(img(:))';
    
 %%
    d=i;
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
 %%
 longitude=54.3569;
    delta_UTC=3.5;
    d=i;
    LT=12;
    LSTM=15*delta_UTC;
    DC=2*pi/365;
    B=DC*(d+10)+0.033*sin(DC*(d-2));
    E0T=(9.87*(sin(2*B))+(7.6)*(sin(B-0.2)));
    TC=4*(longitude-LSTM)+E0T;
    LST=LT+(TC/60);
    tm=LST+1; %tm
    tm=tm- sunrise;
    u=(pi/w)*(ts-tm);
     k=(w/pi)*atan((pi/w)*(ts-tm))
%% calculation mean and amplitude 
 data2=cat(1,lst10,lst13);%Image Read
for i= 1:115116
    y=[data2(:,i)];
     f=[1, cos(((pi/w)*((10- sunrise)-tm))); 1, cos(((pi/w)*((13- sunrise)-tm)))];; 
    a=(inv(transpose(f)*f))*transpose(f)*double(y);
    result(:,i)=a;
end
mean2=result(1,:);
amplitude2=result(2,:);
%%
for i=1:ts
    t=i;
    T =mean+amplitude.*cos((pi/w).*((t)-tm));
   A=amplitude2;
   Thermal_inertia=inertia;
   n=0.5;
   M=ts;
   t=i*3600;
  ww=(2*pi)/(12*3600);
   phase=0;
   %% Calculation of soil heat flux
   G=Thermal_inertia.*(A.*sqrt(n*ww)*sin((n*ww*t)+phase-(pi/8)));
   %% Save the output
    G=reshape(G,318,362); 
   nnaammee=['G_' num2str(ddatte) '_' num2str(i)];
   filenam = ['F:\Thesis\5-G\modis\DTG\G_INERTIA\G\' nnaammee '.tif']
   geotiffwrite( filenam,G,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
end
end