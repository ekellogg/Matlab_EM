calim = zeros(7676,7420,3);
calim(:,:,1) = ReadMRC('../gold_calibration/images/calibration_images/SR_focus.mrc');
calim(:,:,2) = ReadMRC('../gold_calibration/images/calibration_images/SR_1p75.mrc');
calim(:,:,3) = ReadMRC('../gold_calibration/images/calibration_images/SR_3.mrc');

