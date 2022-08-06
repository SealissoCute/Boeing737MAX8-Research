clc
clear
tao = 1.4;
count = 1;
for M = 0:0.1:1
    pe = (0.25*M^2 + (1/24)*M^4*(-tao + 2))*100 ;% change to percent
    pe_all(:,count) =pe;
    count = count+1;
end
plot(0:0.1:1,pe_all)
title('Change in dynamic pressure due to compressibility')
xlabel('Mach number')
ylabel('Percent error in dynamic pressure')