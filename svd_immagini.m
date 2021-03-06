clear all;
close all;

col=gray;
PP = double( imread('Bandiera_italia.jpg') );
[m,n] =size(PP);
P=PP(:,:,1);
[U,S,V] = svd(P); s = diag(S);
k = 1; S(k+1:end,k+1:end) = 0;
PC = U*S*V';
figure(1)
subplot(1,2,1); image(P); title('original image');
colormap(col);
subplot(1,2,2); image(PC); title('compressed (40 S.V.)');
colormap(col);
figure(2)
subplot(1,2,1); image(abs(P-PC)); title('difference');
colormap(col);
subplot(1,2,2); cla; semilogy(s); hold on; plot(k,s(k),'ro');