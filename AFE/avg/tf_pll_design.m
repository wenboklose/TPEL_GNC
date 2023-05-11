kp=49*1;
ki=12*pi*49*2;
tf_pll=tf([kp ki],[1 0])
figure
bode(tf_pll)
grid on