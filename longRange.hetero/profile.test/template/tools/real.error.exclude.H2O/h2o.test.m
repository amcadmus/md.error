beta=2.5;

r01 = r0 - r1;
r02 = r0-r2;
rr01 = norm(r01);
rr02 = norm(r02);
c0 =  0.834;
c1 = -0.417;
c2 = -0.417;

f01r = c0 * c1 * (2 * beta / sqrt(pi) * exp(-beta*beta * rr01 * rr01) + erfc(beta * rr01) / rr01) * r01 / (rr01 * rr01);
f02r = c0 * c2 * (2 * beta / sqrt(pi) * exp(-beta*beta * rr02 * rr02) + erfc(beta * rr02) / rr02) * r02 / (rr02 * rr02);

f01 = c0 * c1 * r01 / (rr01 * rr01 * rr01);
f02 = c0 * c1 * r02 / (rr02 * rr02 * rr02);

f01m = f01 - f01r;
f02m = f02 - f02r;

