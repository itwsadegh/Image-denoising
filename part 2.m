%  part two combined noise
clc 
clear all
 close all
img1=imread('lenna.tif');
% gsimg1=im2double((img1));
 gsimg1=im2double(rgb2gray(img1));
 imshow(gsimg1); title('original')
v=.005; % gaussian noise variance
d=0.005; % salt & pepper noise density 
vs=.005;  % speckle noise variance 
nimg=imnoise(gsimg1,'gaussian',0,v);  % 
nimg1=imnoise(nimg,'salt & pepper',d); %gaussian + salt & pepper 
nimg2=imnoise(nimg,'speckle',vs);      % gaussian + speckle 
nimg3=imnoise(nimg2,'salt & pepper',d);  % gaussian + sant & pepper + speckle 
psr1=psnr(nimg1,gsimg1)
psr2=psnr(nimg2,gsimg1)
psr3=psnr(nimg3,gsimg1)
figure
subplot(1,3,1)
imshow(nimg1);
subplot(1,3,2)
imshow(nimg2);
subplot(1,3,3)
imshow(nimg3);
% wavelet decomposition 
[a1,h1,v1,d1]=dwt2(nimg3,'sym5');
[a2,h2,v2,d2]=dwt2(a1,'sym5');
s1=median(median(abs(h1)))/.6745;
s2=median(median(abs(h2)))/.6745;
[m,n]=size(gsimg1);
dl1=1;
dl2=2;
t1=(s1*sqrt(2*log(m*n)))/(2*dl1-1);
t2=(s2*sqrt(2*log(m*n)))/(2*dl2-1);
hindx1=abs(h1)>=t1;
dindx1=abs(d1)>=t1;
vindx1=abs(v1)>=t1;
hindx2=abs(h2)>=t2;
dindx2=abs(d2)>=t2;
vindx2=abs(v2)>=t2;

%% combined denoising method  
R=0.2;
% level 1
h1hati1=hindx1.*(abs(h1)-(1-R)*t1).*sign(h1);   h1hati2=R*(not(hindx1)).*((abs(h1).^2)/t1);
h1hati=h1hati1+h1hati2;
v1hati1=vindx1.*(abs(v1)-(1-R)*t1).*sign(v1);   v1hati2=R*(not(vindx1)).*((abs(v1).^2)/t1);
v1hati=v1hati1+v1hati2;
d1hati1=dindx1.*(abs(d1)-(1-R)*t1).*sign(d1);    d1hati2=R*(not(dindx1)).*((abs(d1).^2)/t1);
d1hati=d1hati1+d1hati2;
% level 2 
h2hati1=hindx2.*(abs(h2)-(1-R)*t2).*sign(h2);   h2hati2=R*(not(hindx2)).*((abs(h2).^2)/t2);
h2hati=h2hati1+h2hati2;
v2hati1=vindx2.*(abs(v2)-(1-R)*t2).*sign(v2);   v2hati2=R*(not(vindx2)).*((abs(v2).^2)/t2);
v2hati=v2hati1+v2hati2;
d2hati1=dindx2.*(abs(d2)-(1-R)*t2).*sign(d2);    d2hati2=R*(not(dindx2)).*((abs(d2).^2)/t2);
d2hati=d2hati1+d2hati2;
% reconstruction 
a1ni=idwt2(a2,h2hati,v2hati,d2hati,'sym5');
imghati=idwt2(a1ni,h1hati,v1hati,d1hati,'sym5');
% median filter
xi1=medfilt2(imghati);
xi1=medfilt2(xi1);
 peaksnri1=psnr(xi1,gsimg1)
figure 
imshow(xi1); title(' combined denoising output')
%%
xm1=medfilt2(medfilt2(nimg3));
peaksnrm1=psnr(xm1,gsimg1)
%%
figure
subplot(1,2,1)
imshow(nimg3); title('noisy')
subplot(1,2,2)
imshow(imghati) ; title('combined denoising output')

 