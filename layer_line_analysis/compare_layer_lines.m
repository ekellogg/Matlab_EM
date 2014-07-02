%datasets
data = {
'film/tax/my_classums1_69_2400.img', ...
%'film/dyn/my_classums1_69_2400.img', ...
'k2/tax/my_classums1_69_2300.img', ...
%'k2/dyn/my_classums1_69_XX.img'; , ...
};

apix = [1.74 1.74 1.32 1.32];
padfactor = 2;

peaklocs = {};
reso_x = {};
layer_lines = {};

for(i = 1:length(data))
   imgs = ReadImagic(data{i});
   [pks,x,ll] = measure_80A_layer_lines(imgs,apix(i),padfactor);
   peaklocs{i} = pks;
   reso_x{i} = x;
   layer_lines{i} = ll;
end