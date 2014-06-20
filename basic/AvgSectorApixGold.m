gold_line = 2.3 % gold lines show up at 2.3 Angstrom
im_size = 7420;
im_center = 3710; %7420/2
gold_band_ind_begin = 5650;
gold_band_ind_end = 5950;
numimg = 4;

image_index = 1;
sector_size = 30;
ind = 1;
cmap = colormap(jet(length(0:sector_size:360)));
figure()
hold on
for(im = 1:numimg)
sector_apix = [];
subplot(numimg,2,ind)
hold on
cindex=1;
for(th = 0:sector_size:360)
thstart = th*pi/180;
thend = (th+sector_size)*pi/180;
[fr,pw1d,r_ind] = oneDpowerSpectrumAvgSectors(calim(1:7420,:,im),1,thstart,thend);
smo = smoothn(pw1d,35);

search_band_ind_begin = find(r_ind == gold_band_ind_begin);
search_band_ind_end = find(r_ind == gold_band_ind_end);
search_reg = smo(search_band_ind_begin:search_band_ind_end);
semilogy(r_ind(search_band_ind_begin:search_band_ind_end),smo(search_band_ind_begin:search_band_ind_end),...
    '-','Color',cmap(cindex,:),'LineWidth',2)
cindex=cindex+1;
max_ind = r_ind( search_band_ind_begin + find(search_reg == max(search_reg)) - 1 );
apix = 1/((1/gold_line)*(im_size)/(max_ind-im_center))*2
sector_apix = [sector_apix apix];
end
set(gca,'FontSize',16)
set(gca,'FontName','Helvetica')
ind=ind+1;
subplot(numimg,2,ind)
plot(0:sector_size:360,sector_apix,'LineWidth',2)
set(gca,'FontSize',16)
set(gca,'FontName','Helvetica')
ind=ind+1;
end