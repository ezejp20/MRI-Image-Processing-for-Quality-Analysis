%%**************************************************************************
% USE
% Script that produces an ADC map from b value dMRI images and 
% 
% then performs SNR and unifomrity measurements on ADC image
%
% Created by M. Mills and J.Smith 
% Version 3 (06/10/2017)
%
%**************************************************************************

%%

Scanner_details % obtain the scanner details which are relevent

prompt5 = 'What was the room temperature? Enter in degrees C';
    Temp = input(prompt5);

end_excel_ap % close current workbooks so data can be written in

%%

% warning('off','curvefit:fit:nonDoubleYData');
warning('off','all');

 
 fprintf('Select the b_50 image');
 [file_name, path_name, filter_index] = uigetfile('*.*');
 b1 = dicomread(strcat(path_name,file_name));

 fprintf('Select the b_400 image');
 [file_name, path_name, filter_index] = uigetfile('*.*');
 b2 = dicomread(strcat(path_name,file_name));
 
 fprintf('Select the b_800 image');
 [file_name, path_name, filter_index] = uigetfile('*.*');
 b3 = dicomread(strcat(path_name,file_name));
 
 fprintf('Select the b_1000 image');
 [file_name, path_name, filter_index] = uigetfile('*.*');
 b4 = dicomread(strcat(path_name,file_name));
 
 % info = dicominfo(b1);
 ADC = zeros(size(b1));

ADC = b1 ;
t = size(b1);
tt = t(1)
goftab = zeros(size(b1));
x = [50;400;800;1000];

Coval = zeros(size(b1));

for i=1:tt
    for j = 1:tt
      
        y = [b1(i,j);b2(i,j);b3(i,j);b4(i,j)];
        
        [f,gof] = fit(x,y,'exp1');
        
        coval=coeffvalues(f);
        
        Coval(i,j) = coval(2);
        
        %goftab(i,j) = gof;
    
        ADC(i,j) = abs(Coval(i,j)*1000000);
        
        
    end
    i

end

figure,imshow(ADC, 'DisplayRange', [min(min(ADC)) max(mean(ADC))]);

disp('draw a line from the top of the phantom to the bottom of the phantom (longest line) :\n');


hy = imline; % creates line
posy = wait(hy); % holds for interactive placement
posy = getPosition(hy); % gets first and last coord of line [x1 y1; x2 y2];
y1 = posy(1,2); 
y2 = posy(2,2);

if y1 < y2
    yleng = y2 - y1;
    ymid = y2 - (0.5*yleng);
else
    yleng = y1 - y2;
    ymid = y1 - (0.5*yleng);
end

disp('draw a line from the L --> R edge of the phantom :\n');
hx = imline; % creates line
posx = wait(hx); % holds for interactive placement
posx = getPosition(hx); % gets first and last coord of line [x1 y1; x2 y2];
x1 = posx(1,1); 
x2 = posx(2,1);

if x1 < x2
    xleng = x2 - x1;
    xmid = x2 - (0.5*xleng);
else
    xleng = x1 - x2;
    xmid = x1 - (0.5*xleng);
end


if xleng < yleng
    width = 0.86*xleng;
    height = 0.86*xleng;
    xmin = xmid - (0.43*xleng);
    ymin = ymid - (0.43*xleng);
else
    width = 0.86*yleng;
    height = 0.86*yleng;
    xmin = xmid - (0.43*yleng);
    ymin = ymid - (0.43*yleng);
end



posBIG = [xmin ymin width height];

fprintf('Drag ROI, it will snap to location and correct size :\n'); % select roi - do not close image
BIG = imellipse; % creates an ellipse
setPosition(BIG, posBIG); % fits ellipse to coords

widthSm = 0.2*width; % ROIs inside are 20% diameter of large ROI
heightSm = 0.2*height;

possmall = [xmin ymin widthSm heightSm];


 pos1 = possmall; % creates new vector for expanded ROI
    pos2 = possmall;
    pos3 = possmall;
    pos4 = possmall;
    pos5 = possmall;

    pos1(1) = 0.5*posBIG(3) - 0.5*possmall(3) + posBIG(1);
    pos1(2) = 0.2*posBIG(4) - 0.5*possmall(4) + posBIG(2);
    pos2(1) = 0.2*posBIG(3) - 0.5*possmall(3)+ posBIG(1);
    pos2(2) = 0.5*posBIG(4) - 0.5*possmall(4)+ posBIG(2);
    pos3(1) = pos1(1);
    pos3(2) = pos2(2);
    pos4(1) = 0.8*posBIG(3) - 0.5*possmall(3)+ posBIG(1);
    pos4(2) = pos2(2);
    pos5(1) = pos1(1);
    pos5(2) = 0.8*posBIG(4) - 0.5*possmall(4) + posBIG(2);

    fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
    h1 = imellipse;
    setPosition(h1, pos1);

    fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
    h2 = imellipse;
    setPosition(h2, pos2);

    fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
    h3 = imellipse;
    setPosition(h3, pos3);

    fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
    h4 = imellipse;
    setPosition(h4, pos4);

    fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
    h5 = imellipse;
    setPosition(h5, pos5);

    BW1 = createMask(h1);
    BW2 = createMask(h2);
    BW3 = createMask(h3); 
    BW4 = createMask(h4);
    BW5 = createMask(h5);
 
    A1=ADC(BW1);
    A2=ADC(BW2);
    A3=ADC(BW3);
    A4=ADC(BW4);
    A5=ADC(BW5);
    C1=[];
    C2=[];
    C3=[];
    C4=[];
    C5=[];
for i=1:length(A1)
   C1(i)=A1(i); 
end
for i=1:length(A2)
   C2(i)=A2(i); 
end
for i=1:length(A3)
   C3(i)=A3(i); 
end
for i=1:length(A4)
   C4(i)=A4(i); 
end
for i=1:length(A5)
   C5(i)=A5(i); 
end
 
    mean_signal_1 = mean(C1);
    mean_signal_2 = mean(C2);
    mean_signal_3 = mean(C3);
    mean_signal_4 = mean(C4);
    mean_signal_5 = mean(C5); 

    ADC_Av_mean_signal = (mean_signal_1+mean_signal_2+mean_signal_3+mean_signal_4+mean_signal_5)/5;
    
    stdev_1 = std(C1);
    stdev_2 = std(C2);
    stdev_3 = std(C3);
    stdev_4 = std(C4);
    stdev_5 = std(C5);
    Av_stdev_signal = (stdev_1+stdev_2+stdev_3+stdev_4+stdev_5)/5;
 
Va = xleng*0.4;
Va = round(Va);
        
    figure
h = imshow(ADC, 'DisplayRange', [min(min(ADC)) max(mean(ADC))]);

fprintf('Draw Line along frequency encoding direction (L->R)'); % select roi - do not close image
[x,y] = ginput(2);
xline = [x(1),x(2)];
yline = [y(1), y(1)];
a = improfile(ADC, xline, yline);
fe1 = a( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Repeat for FE');
[x2,y2] = ginput(2);
xline2 = [x2(1),x2(2)];
yline2 = [y2(1), y2(1)];
b = improfile(ADC, xline2, yline2); % draw a horizontal line across image
fe2 = b( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Repeat for FE');
[x3,y3] = ginput(2);
xline3 = [x3(1),x3(2)];
yline3 = [y3(1), y3(1)];
c = improfile(ADC, xline3, yline3); % draw a horizontal line across image
fe3 = c( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Repeat for FE');
[x4,y4] = ginput(2);
xline4 = [x4(1),x4(2)];
yline4 = [y4(1), y4(1)];
d = improfile(ADC, xline4, yline4); % draw a horizontal line across image
fe4 = d( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Final FE');
[x5,y5] = ginput(2);
xline5 = [x5(1),x5(2)];
yline5 = [y5(1), y5(1)];
e = improfile(ADC, xline5, yline5); % draw a horizontal line across image
fe5 = e( (ceil(end/2))-Va : (ceil(end/2))+Va);

FE = [fe1 fe2 fe3 fe4 fe5];
FEav = mean(FE,2);
Va = yleng*0.4;
Va = round(Va);

fprintf('\n Draw Line along phase encoding direction (top->bottom)'); % select roi - do not close image
[x6,y6] = ginput(2);
xline6 = [x6(1),x6(1)];
yline6 = [y6(1), y6(2)];
f = improfile(ADC, xline6, yline6); % draw a horizontal line across image
pe1 = f( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Repeat for PE'); % select roi - do not close image
[x7,y7] = ginput(2);
xline7 = [x7(1),x7(1)];
yline7 = [y7(1), y7(2)];
g = improfile(ADC, xline7, yline7); % draw a horizontal line across image
pe2 = g( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Repeat for PE'); % select roi - do not close image
[x8,y8] = ginput(2);
xline8 = [x8(1),x8(1)];
yline8 = [y8(1), y8(2)];
h = improfile(ADC, xline8, yline8); % draw a horizontal line across image
pe3 = h( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Repeat for PE'); % select roi - do not close image
[x9,y9] = ginput(2);
xline9 = [x9(1),x9(1)];
yline9 = [y9(1),y9(2)];
j = improfile(ADC, xline9, yline9); % draw a horizontal line across image
pe4 = j( (ceil(end/2))-Va : (ceil(end/2))+Va);

fprintf('\n Final PE'); % select roi - do not close image
[x10,y10] = ginput(2);
xline10 = [x10(1),x10(1)];
yline10 = [y10(1), y10(2)];
k = improfile(ADC, xline10, yline10); % draw a horizontal line across image
pe5 = k( (ceil(end/2))-Va : (ceil(end/2))+Va);

PE = [pe1 pe2 pe3 pe4 pe5];
PEav = mean(PE,2);

FEmode = mode(FEav);
FElower = FEmode - (0.1*FEmode);
FEupper = FEmode + (0.1*FEmode);
FEop = FEav ( FEav>FElower & FEav<FEupper);
FractionalUniformityFE = 100* (size(FEop))/size(FEav);

PEmode = mode(PEav);
PElower = PEmode - (0.1*PEmode);
PEupper = PEmode + (0.1*PEmode);
PEop = PEav ( PEav>PElower & PEav<PEupper);
FractionalUniformityPE = 100* (size(PEop))/size(PEav);

ADC_FractionalUniformityAv = (FractionalUniformityFE + FractionalUniformityPE)/2;


ADC_Av_mean_signal ; 

ADC_FractionalUniformityAv ; 

%% now SNR and UNIV

%[I, path] = uigetfile('*.*');
%info = dicominfo(fullfile(path,I));
%I = dicomread(info);
%I = imread(I);

% if not dicom image substitute previous with I = imread('image.jpg');
% [I2, path2] = uigetfile('*.*');
% info = dicominfo(fullfile(path2,I2));
% I2 = dicomread(info);
% %I2 = imread(I2);
% 
% figure
% h = imshow(I, 'DisplayRange', [min(min(I)) max(max(I))]);
% 
% disp('draw a line from the top of the phantom to the bottom of the phantom (longest line) :\n');
% 
% figure
% h = imshow(I, 'DisplayRange', [min(min(I)) max(max(I))]);
% 
% hy = imline; % creates line
% posy = wait(hy); % holds for interactive placement
% posy = getPosition(hy); % gets first and last coord of line [x1 y1; x2 y2];
% y1 = posy(1,2); 
% y2 = posy(2,2);
% 
% if y1 < y2
%     yleng = y2 - y1;
%     ymid = y2 - (0.5*yleng);
% else
%     yleng = y1 - y2;
%     ymid = y1 - (0.5*yleng);
% end
% 
% disp('draw a line from the L --> R edge of the phantom :\n');
% hx = imline; % creates line
% posx = wait(hx); % holds for interactive placement
% posx = getPosition(hx); % gets first and last coord of line [x1 y1; x2 y2];
% x1 = posx(1,1); 
% x2 = posx(2,1);
% 
% if x1 < x2
%     xleng = x2 - x1;
%     xmid = x2 - (0.5*xleng);
% else
%     xleng = x1 - x2;
%     xmid = x1 - (0.5*xleng);
% end
% 
% 
% 
% if xleng < yleng
%     width = 0.86*xleng;
%     height = 0.86*xleng;
%     xmin = xmid - (0.43*xleng);
%     ymin = ymid - (0.43*xleng);
% else
%     width = 0.86*yleng;
%     height = 0.86*yleng;
%     xmin = xmid - (0.43*yleng);
%     ymin = ymid - (0.43*yleng);
% end
% 
% 
% 
% posBIG = [xmin ymin width height];
% 
% fprintf('Drag ROI, it will snap to location and correct size :\n'); % select roi - do not close image
% BIG = imellipse; % creates an ellipse
% setPosition(BIG, posBIG); % fits ellipse to coords
% 
% widthSm = 0.2*width; % ROIs inside are 20% diameter of large ROI
% heightSm = 0.2*height;
% 
% possmall = [xmin ymin widthSm heightSm];
% 
% 
%  pos1 = possmall; % creates new vector for expanded ROI
%     pos2 = possmall;
%     pos3 = possmall;
%     pos4 = possmall;
%     pos5 = possmall;
% 
%     pos1(1) = 0.5*posBIG(3) - 0.5*possmall(3) + posBIG(1);
%     pos1(2) = 0.2*posBIG(4) - 0.5*possmall(4) + posBIG(2);
%     pos2(1) = 0.2*posBIG(3) - 0.5*possmall(3)+ posBIG(1);
%     pos2(2) = 0.5*posBIG(4) - 0.5*possmall(4)+ posBIG(2);
%     pos3(1) = pos1(1);
%     pos3(2) = pos2(2);
%     pos4(1) = 0.8*posBIG(3) - 0.5*possmall(3)+ posBIG(1);
%     pos4(2) = pos2(2);
%     pos5(1) = pos1(1);
%     pos5(2) = 0.8*posBIG(4) - 0.5*possmall(4) + posBIG(2);
% 
%     fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
%     h1 = imellipse;
%     setPosition(h1, pos1);
% 
%     fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
%     h2 = imellipse;
%     setPosition(h2, pos2);
% 
%     fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
%     h3 = imellipse;
%     setPosition(h3, pos3);
% 
%     fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
%     h4 = imellipse;
%     setPosition(h4, pos4);
% 
%     fprintf('\n Draw any ROI - ROI will autoplace'); % select roi - do not close image
%     h5 = imellipse;
%     setPosition(h5, pos5);
% 
%     BW1 = createMask(h1);
%     BW2 = createMask(h2);
%     BW3 = createMask(h3); 
%     BW4 = createMask(h4);
%     BW5 = createMask(h5);
%  
%     A1=I(BW1);
%     A2=I(BW2);
%     A3=I(BW3);
%     A4=I(BW4);
%     A5=I(BW5);
%     C1=[];
%     C2=[];
%     C3=[];
%     C4=[];
%     C5=[];
% for i=1:length(A1)
%    C1(i)=A1(i); 
% end
% for i=1:length(A2)
%    C2(i)=A2(i); 
% end
% for i=1:length(A3)
%    C3(i)=A3(i); 
% end
% for i=1:length(A4)
%    C4(i)=A4(i); 
% end
% for i=1:length(A5)
%    C5(i)=A5(i); 
% end
%  
%     mean_signal_1 = mean(C1);
%     mean_signal_2 = mean(C2);
%     mean_signal_3 = mean(C3);
%     mean_signal_4 = mean(C4);
%     mean_signal_5 = mean(C5); 
% 
%     
%  
%     stdev_1 = std(C1);
%     stdev_2 = std(C2);
%     stdev_3 = std(C3);
%     stdev_4 = std(C4);
%     stdev_5 = std(C5);
%     Av_stdev_signal = (stdev_1+stdev_2+stdev_3+stdev_4+stdev_5)/5;
%  
%  
%     Inoise = imsubtract(I, I2); % subtracts each element of array from one another
%  
%     figure
%     h_noise = imshow(Inoise, 'DisplayRange', [min(min(Inoise)) max(max(Inoise))]);
%  
%  
%     A1=Inoise(BW1);
%     A2=Inoise(BW2);
%     A3=Inoise(BW3);
%     A4=Inoise(BW4);
%     A5=Inoise(BW5);
%     D1=[];
%     D2=[];
%     D3=[];
%     D4=[];
%     D5=[];
% for i=1:length(A1)
%    D1(i)=A1(i); 
% end
% for i=1:length(A2)
%    D2(i)=A2(i); 
% end
% for i=1:length(A3)
%    D3(i)=A3(i); 
% end
% for i=1:length(A4)
%    D4(i)=A4(i); 
% end
% for i=1:length(A5)
%    D5(i)=A5(i); 
% end
%  
%     stdev_noise_1 = std(D1);
%     stdev_noise_2 = std(D2);
%     stdev_noise_3 = std(D3);
%     stdev_noise_4 = std(D4);
%     stdev_noise_5 = std(D5);
%  
%     SNR1 = ((sqrt(2))*mean_signal_1/stdev_noise_1);
%     SNR2 = ((sqrt(2))*mean_signal_2/stdev_noise_2);
%     SNR3 = ((sqrt(2))*mean_signal_3/stdev_noise_3);
%     SNR4 = ((sqrt(2))*mean_signal_4/stdev_noise_4);
%     SNR5 = ((sqrt(2))*mean_signal_5/stdev_noise_5);
%     SNRav = (SNR1+SNR2+SNR3+SNR4+SNR5)/5;
%     
% Va = xleng*0.4;
% Va = round(Va);
% 
% figure
% h = imshow(I, 'DisplayRange', [min(min(I)) max(max(I))]);
% 
% fprintf('Draw Line along frequency encoding direction (L->R)'); % select roi - do not close image
% [x,y] = ginput(2);
% xline = [x(1),x(2)];
% yline = [y(1), y(1)];
% a = improfile(I, xline, yline);
% fe1 = a( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Repeat for FE');
% [x2,y2] = ginput(2);
% xline2 = [x2(1),x2(2)];
% yline2 = [y2(1), y2(1)];
% b = improfile(I, xline2, yline2); % draw a horizontal line across image
% fe2 = b( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Repeat for FE');
% [x3,y3] = ginput(2);
% xline3 = [x3(1),x3(2)];
% yline3 = [y3(1), y3(1)];
% c = improfile(I, xline3, yline3); % draw a horizontal line across image
% fe3 = c( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Repeat for FE');
% [x4,y4] = ginput(2);
% xline4 = [x4(1),x4(2)];
% yline4 = [y4(1), y4(1)];
% d = improfile(I, xline4, yline4); % draw a horizontal line across image
% fe4 = d( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Final FE');
% [x5,y5] = ginput(2);
% xline5 = [x5(1),x5(2)];
% yline5 = [y5(1), y5(1)];
% e = improfile(I, xline5, yline5); % draw a horizontal line across image
% fe5 = e( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% FE = [fe1 fe2 fe3 fe4 fe5];
% FEav = mean(FE,2);
% Va = yleng*0.4;
% Va = round(Va);
% 
% fprintf('\n Draw Line along phase encoding direction (top->bottom)'); % select roi - do not close image
% [x6,y6] = ginput(2);
% xline6 = [x6(1),x6(1)];
% yline6 = [y6(1), y6(2)];
% f = improfile(I, xline6, yline6); % draw a horizontal line across image
% pe1 = f( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Repeat for PE'); % select roi - do not close image
% [x7,y7] = ginput(2);
% xline7 = [x7(1),x7(1)];
% yline7 = [y7(1), y7(2)];
% g = improfile(I, xline7, yline7); % draw a horizontal line across image
% pe2 = g( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Repeat for PE'); % select roi - do not close image
% [x8,y8] = ginput(2);
% xline8 = [x8(1),x8(1)];
% yline8 = [y8(1), y8(2)];
% h = improfile(I, xline8, yline8); % draw a horizontal line across image
% pe3 = h( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Repeat for PE'); % select roi - do not close image
% [x9,y9] = ginput(2);
% xline9 = [x9(1),x9(1)];
% yline9 = [y9(1),y9(2)];
% j = improfile(I, xline9, yline9); % draw a horizontal line across image
% pe4 = j( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% fprintf('\n Final PE'); % select roi - do not close image
% [x10,y10] = ginput(2);
% xline10 = [x10(1),x10(1)];
% yline10 = [y10(1), y10(2)];
% k = improfile(I, xline10, yline10); % draw a horizontal line across image
% pe5 = k( (ceil(end/2))-Va : (ceil(end/2))+Va);
% 
% PE = [pe1 pe2 pe3 pe4 pe5];
% PEav = mean(PE,2);
% 
% FEmode = mode(FEav);
% FElower = FEmode - (0.1*FEmode)
% FEupper = FEmode + (0.1*FEmode);
% FEop = FEav ( FEav>FElower & FEav<FEupper);
% FractionalUniformityFE = 100* (size(FEop))/size(FEav);
% 
% PEmode = mode(PEav);
% PElower = PEmode - (0.1*PEmode)
% PEupper = PEmode + (0.1*PEmode);
% PEop = PEav ( PEav>PElower & PEav<PEupper);
% FractionalUniformityPE = 100* (size(PEop))/size(PEav);
% 
% FractionalUniformityAv = (FractionalUniformityFE + FractionalUniformityPE)/2;


%%
 sheet = 'Routine_QA_pros';
 QA = xlsread('MRQA_MATLAB.xlsx', sheet, 'A1:Z500');
 t = {Date Scanner Tes Temp ADC_Av_mean_signal ADC_FractionalUniformityAv} ; % t = {Date Scanner Tes Temp ADC_Av_mean_signal ADC_FractionalUniformityAv SNRav FractionalUniformityAv} ;
 nRows = (size(QA,1));
 nRows = nRows + 3;
 K = num2str(nRows);
 c = strcat('A',K);
 xlswrite('MRQA_MATLAB',t, sheet,c);






 
  
