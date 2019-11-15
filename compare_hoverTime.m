close all; clear;

name = categorical({'GTQ Mini' 'IRIS+ (with payload)' 'PD6-AW2 BASIC' 'Fox 4'});
name = reordercats(name,{'GTQ Mini' 'IRIS+ (with payload)' 'PD6-AW2 BASIC' 'Fox 4'});
manufacturer = [6.8 22 30 25]';
mine = [4.9 21.1 35.6 27.8]';

display(round(((manufacturer-mine)./manufacturer)*100)/100);

figure;
bar(name, [mine manufacturer]);
title('Comparison of estimated hover time');
legend('Methodology result', 'Reference value', 'Orientation', 'horizontal', 'Location', 'northoutside');
ylabel('Time [min]'); grid;