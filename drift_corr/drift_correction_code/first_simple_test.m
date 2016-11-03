setup_env
dfcorr = simple_drift_correction(fr);
dfsum = cumsum(dfcorr,1);
subplot(3,2,1); plot(dfsum(:,1),dfsum(:,2));

dfsum = gather(dfsum);
summed_im = apply_shifts_to_stack_v2(dfsum,fr);
subplot(3,2,3);
showImage(medfilt2(summed_im, [5 5]));
subplot(3,2,4);
showImage(medfilt2(sum(fr,3),[5 5]))
subplot(3,2,5);
tt = (medfilt2(fftim(summed_im),[5 5]));
showImage(tt);
subplot(3,2,6)
tt = (medfilt2(fftim(sum(fr,3)),[5 5]));
showImage(tt);
subplot(3,2,2)
[testf,testavg,test_ind] = oneDpowerSpectrum(summed_im,1.32);
[origf,origavg,orig_ind] = oneDpowerSpectrum(sum(fr,3),1.32);
semilogy(testf,testavg)
hold on
semilogy(origf,origavg,'-red')
legend('simple-drift-corr','sum')