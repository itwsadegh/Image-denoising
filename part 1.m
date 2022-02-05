% part 1 gaussian noise 
clc 
clear all
 close all
img1=imread('lenna.tif');
gsimg1=im2double(rgb2gray(img1));
% gsimg1=imresize(gsimg1,0.5);
% imshow(img1); title('lenna')
V=.00151;    % noise variance 
nimg=imnoise(gsimg1,'gaussian',0,V); % noisy image 
psr=psnr(nimg,gsimg1)
figure;
subplot(1,2,1); imshow(gsimg1);
 subplot(1,2,2);   imshow(nimg);
% wavelet decomposition 
[a1,h1,v1,d1]=dwt2(nimg,'sym5');
[a2,h2,v2,d2]=dwt2(a1,'sym5');
s1=median(median(abs(h1)))/.6745;
 s2=median(median(abs(h2)))/.6745;
[m,n]=size(gsimg1);
% [p,q]=size(h2);
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
%% hard threshold
% level 1
 h1hat=hindx1.*(h1);
 d1hat=dindx1.*(d1);
 v1hat=vindx1.*(v1);
% level 2
h2hat=hindx2.*(h2);
 d2hat=dindx2.*(d2);
 v2hat=vindx2.*(v2);
% reconstruction 
a1n=idwt2(a2,h2hat,v2hat,d2hat,'sym5');
imghat=idwt2(a1n,h1hat,v1hat,d1hat,'sym5');
%
figure;
imshow(imghat,[]); title('hard  threshold output')
peaksnr=psnr(imghat,gsimg1)
%% soft threshold 
% level 1
 h1hats=sign(h1).*hindx1.*(abs(h1)-t1);
 d1hats=sign(d1).*dindx1.*(abs(d1)-t1);
 v1hats=sign(v1).*vindx1.*(abs(v1)-t1);
% level 2
 h2hats=sign(h2).*hindx2.*(abs(h2)-t2);
 d2hats=sign(d2).*dindx2.*(abs(d2)-t2);
 v2hats=sign(v2).*vindx2.*(abs(v2)-t2);
% reconstruction 
a1ns=idwt2(a2,h2hats,v2hats,d2hats,'sym5');
imghats=idwt2(a1ns,h1hats,v1hats,d1hats,'sym5');
figure 
imshow(imghats);title('soft threshold output')
peaksnrs=psnr(imghats,gsimg1)
%% improved threshold 
R=0.1;
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
figure 
imshow(imghati); title('improved threshold output')
peaksnri=psnr(imghati,gsimg1)
%% median filter 
imgm1=medfilt2(medfilt2(nimg));
psnrm=psnr(imgm1,gsimg1)
imshow(imgm1) ; title('median output')
%%
subplot(1,6,1)
imshow(gsimg1)
title('original')
subplot(1,6,2)
imshow(nimg)
title('noisy')
subplot(1,6,3)
imshow(imghat)
title('hard')
subplot(1,6,4)
imshow(imghats)
title('soft')
subplot(1,6,5)
imshow(imghati)
title('improved')
subplot(1,6,6)
imshow(imgm1)
title('median')


