function out = RK(yn)
global h;

kk1 = wholesystem(yn);
kk2 = wholesystem(yn + kk1*h/2);
kk3 = wholesystem(yn + kk2*h/2);
kk4 = wholesystem(yn + kk3*h);

out = yn + (kk1 + 2*kk2 + 2*kk3 + kk4)*h/6;

end