%% Calculation of Soil Heat Flux at Night using MODIS images
% In this code, the energy balance equation is used to estimate the soil heat flux.

%!! Important Note!!
% This code is prepared to estimate soil heat flux at night. To use it, time series images of land surface temperature at night, emissivity, ndvi vegetation cover index and albedo are needed.
% MODIS sensor has these images as products. Before running the code, it is necessary to resize the images.
% Pay attention to the naming and numbering order of the images inside the folders.

%%!!Code execution steps
%1-Enter the path to save LST, Albedo, Gasilmandi and NDVI images
%2-Enter the number of images you want in front of for.
%3-Change the output path and save the image

% This code was written by Fahime Arabi. if you have any questions about it, I will answer you with the following email:
%Fahimearabi1993@gmail.com
%%
clc;
clear;
close all;
%Read Images from my computer
path_lst='F:\Thesis\5-G\modis\lst\2020\MYD\NIGHT-RESIZE'
path_ndvi='F:\Thesis\5-G\modis\ndvi\5-ndvi-daily\2020'
path_albedo='F:\Thesis\5-G\modis\albedo\image\albedo2020'
path_emi='F:\Thesis\5-G\modis\emissivity\2020'
%%
for i=1:365  %Enter the number of images you want to process
    cd(path_ndvi);
    a=dir('*.tif');
    namendvi=a(i).name;
    [img,Ref]=geotiffread([path_ndvi, '\', namendvi]);
    info=geotiffinfo([path_ndvi, '\', namendvi]);
    ndvi=img(:,:,:);

    cd(path_lst);
    a=dir('*.tif');
    namelst=a(i).name;
    [img,Ref]=geotiffread([path_lst, '\', namelst]);
     [B4,RefMatrx]=geotiffread(namelst);
    InfoB4=geotiffinfo('NDVI_doy2014001.tif');
    InfoB4=geotiffinfo([path_lst, '\', namelst]);
    lst=img(:,:,:);
    
    cd(path_emi);
    a=dir('*.tif');
    nameemi=a(i).name;
    [img,Ref]=geotiffread([path_emi, '\', nameemi]);
    info=geotiffinfo([path_emi, '\', nameemi]);
    emi=img(:,:,:);
    
    cd(path_albedo);
    a=dir('*.tif');
    namealbedo=a(i).name;
    [img,Ref]=geotiffread([path_albedo, '\', namealbedo]);
    strname=strsplit(namealbedo,'_');
    ddatte=strname{1,2};
    info=geotiffinfo([path_albedo, '\', namealbedo]);
    albedo=img(:,:,:);
   disp('Images read');
    %% LAI index
  lai=0.57 .*(exp(2.33 .* ndvi)); %Leaf area index (LAI)
%% emissivity εNB & ε0
  kk=size(img);
  difx=kk(1,2); 
  dify=kk(1,1);
  images2=double(zeros(dify,difx));
  ep_0=images2*0.98;
  c=zeros(dify,difx);
  c(ndvi>0.1 & ndvi<0.6)=1;
      ss1 =((0.971 .*(1-((ndvi-0.1)./(0.5))))+(0.987 .*((ndvi-0.1) ./(0.5)))) .* c;
      tt1 =((0.977 .*(1-((ndvi-0.1)./(0.5))))+(0.989 .*((ndvi-0.1) ./(0.5)))) .* c;
      
c=zeros(dify,difx);
c(ndvi<=0.1)=1;
     ss2= 0.971 .* c;
     tt2= 0.977 .* c;

c=zeros(dify,difx);
c(ndvi>=0.6)=1;
     ss3=0.987 .* c;
     tt3=0.989 .* c;
     
c=zeros(dify,difx);
c(ndvi>0 & lai<3)=1;
y1=0.95+((0.01).*ndvi);
y1=y1 .* c;

c=zeros(dify,difx);
c(ndvi>0 & lai>=3)=1;
y2=0.98 .* c;

c=zeros(dify,difx);
c(ndvi<=0)=1;
y3=0.985 .* c;

ep_0=y1+y2+y3;

all_pixel=difx*dify;
lst_ch=(lst(:))';
ndv_ch=(ndvi(:))';
lst_sort=sort(lst_ch);
ndv_sort=sort(ndv_ch);

images2=double(zeros(dify,difx));
ep_0=emi;
c=zeros(dify,difx);
c(ndvi>0.1 & ndvi<0.6)=1;
      ss1 =((0.971 .*(1-((ndvi-0.1)./(0.5))))+(0.987 .*((ndvi-0.1) ./(0.5)))) .* c;
      tt1 =((0.977 .*(1-((ndvi-0.1)./(0.5))))+(0.989 .*((ndvi-0.1) ./(0.5)))) .* c;
      
c=zeros(dify,difx);
c(ndvi<=0.1)=1;
     ss2= 0.971 .* c;
     tt2= 0.977 .* c;

c=zeros(dify,difx);
c(ndvi>=0.6)=1;
     ss3=0.987 .* c;
     tt3=0.989 .* c;
     
c=zeros(dify,difx);
c(ndvi>0 & lai<3)=1;
y1=0.95+((0.01).*ndvi);
y1=y1 .* c;

c=zeros(dify,difx);
c(ndvi>0 & lai>=3)=1;
y2=0.98 .* c;

c=zeros(dify,difx);
c(ndvi<=0)=1;
y3=0.985 .* c;

ep_0=y1+y2+y3;
%%  Rn_net radiation flux
rl_out=(ep_0).*((5.68).*10^-8).*((lst).^4);
rs_in=0;
rl_in=0;
Rn=((1-albedo).*(rs_in))+(rl_in)-(rl_out)-((1-ep_0).*(rl_in));
%% Soil Heat Flux (G)
G=(((lst-273.15)./albedo).*((((albedo).^2)*0.0074)+((albedo).*(0.0038))).*(1-((ndvi.^4).*(0.98)))).* Rn;
%% Save the output
nnaammee=['Gn' num2str(ddatte)];
    filenam = ['F:\Thesis\5-G\modis\soil heat flux\2020\G-MYD\NIGHT\' nnaammee '.tif']
geotiffwrite( filenam,G,RefMatrx,'GeoKeyDirectoryTag',InfoB4.GeoTIFFTags.GeoKeyDirectoryTag);
end