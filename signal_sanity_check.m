clear
clc

rx_signal_1 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_new_signal_05_16_2020/antenna1.bin',10000000000000000);

rx_signal_2 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_uncalibrated/antenna1.bin',10000000000000000);

rx_signal_3 = read_complex_binary('/Users/changhange/Downloads/data_0521/antenna_calibrated/antenna1.bin',10000000000000000);

rx_signal_1 = rx_signal_1(1:100000);
rx_signal_2 = rx_signal_2(1:100000);
rx_signal_3 = rx_signal_3(1:100000);

rx_signal_1_corr = xcorr(rx_signal_1, rx_signal_1);
rx_signal_2_corr = xcorr(rx_signal_2, rx_signal_2);
rx_signal_3_corr = xcorr(rx_signal_3, rx_signal_3);

rx_signal_1_corr = rx_signal_1_corr(length(rx_signal_1):end);
rx_signal_2_corr = rx_signal_2_corr(length(rx_signal_1):end);
rx_signal_3_corr = rx_signal_3_corr(length(rx_signal_1):end);

subplot(3,1,1);plot(abs(rx_signal_1_corr))
subplot(3,1,2);plot(abs(rx_signal_2_corr))
subplot(3,1,3);plot(abs(rx_signal_3_corr))