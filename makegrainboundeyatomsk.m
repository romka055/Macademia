function l=makegrainboundeyatomsk(ang,a);
%%%
% 26/3/2020
% This function makes a periodic grain boundary for lammps using Atomsk
% "ang" is the angle of the grain boundary
% "a" is the lattice spacing. This version works for 100/010/001 crystals
% only
% The main idea is that you need to find the periodicity in the scewed
% crystal along the X direction and it is based on arithmetics of whole
% numbers.

filename = [num2str(ang),'.sh'];
i=0;
f=50;
b=a;

res = 100;
while abs(res)>0.05
f=f+1;
mn=tand(ang)*a/b;
[m,n]=rat(round(mn*f)/f);
x=tand(ang)*n/sind(ang)*a;
y=n*sind(ang)*a;
%%%
[m1,n1]=rat(mn);
x1=tand(ang)*n1/sind(ang)*a;
y1=n1*sind(ang)*a;
res=round(x1/x)-x1/x;
i=i+1;
if i>100
break;
end
end
xres=x*ceil((20*a)/x);
yres=y*ceil((20*a)/y);
l=[xres yres];

% Now lets write the file
fid=fopen(filename,'w');
fprintf(fid,'#!/bin/bash \nrm -f box.txt\n');
fprintf(fid,'echo "box %f %f %f" > box.txt\n',xres,yres,a);
fprintf(fid,'echo "node 0.5*box 0.25*box 0 0° 0° -%f°" >> box.txt\n',ang);
fprintf(fid,'echo "node 0.5*box 0.75*box 0 0° 0° %f°" >> box.txt\n',ang);
fprintf(fid,'rm -f Ta*.xsf\nrm -f Ta*.dat\nrm -f %d.cfg\n rm -f Ta.cfg\n',ang);
fprintf(fid,'./atomsk --create bcc %f Ta orient [100] [010] [001] Ta.xsf\n./atomsk --polycrystal Ta.xsf box.txt Ta.cfg\n./atomsk Ta.cfg -wrap -remove-doubles 1 %d.cfg\n',a,ang);
fclose(fid);
end
