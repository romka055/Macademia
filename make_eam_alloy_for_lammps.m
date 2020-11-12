%% S script to make EAM/Alloy pair potential for Lammps
% This is an example of Zinc based on Belashchenko - High Temperature, 2012, Vol. 50, No. 1, pp. 61–69.
% U = F(psi(r))+phi(r)
% Note that for lammps (eam/alloy) you need to input phi(r)*r
clear;
clc;
% EAM/Alloy generator.
% Created by Roman Kositski 2013
%%% ----------------------------------
Nrho=5000;
rhomax=10;
drho=rhomax/Nrho;
Nr=5000;
rmax=7.6;
rmin=0;
dr=(rmax-rmin)/(Nr);
rl=linspace(rmin,rmax-dr,Nr);
mass=65.38;
lattice_constant=5.41;
lattice_type='FCC';
Z=30; % Atomic number
qe=1; % Electron Charge
rb=25; % Bohr Radius - Biersack, J. P., and Ziegler, J. F., 1982, Nucl. Instrum. Meth., 141, 93.
%%% ------------------------------
% Start writing to file
filename = 'zinc.eam.alloy';
fid = fopen(filename, 'w');
% ------------------------------
% Alloy file header
%lines 1,2,3 = comments (ignored)
fprintf (fid,'LAMMPS Potential File - Liquid Zinc\n');
fprintf (fid,'Belashchenko - High Temperature, 2012, Vol. 50, No. 1, pp. 61–69. \n');
fprintf (fid,'Implemented by Roman Kositski 2020 \n');
%line 4: Nelements Element1 Element2 ... ElementN
fprintf (fid,'1 Zn\n');
%line 5: Nrho, drho, Nr, dr, cutoff
fprintf (fid, '%d %12.6e %d %12.6e %12.6e\n',Nrho, drho, Nr, dr, rmax);
% ----------------------------------
%%%% Nelements sections, one for each element  %%%%
% line 1 = atomic number, mass, lattice constant, lattice type (e.g. FCC)
fprintf (fid,'%d %12.6f %12.6f %s \n',Z, mass, lattice_constant, lattice_type);
%%======
%embedding function F(rho) (Nrho values)

rhol=linspace(0,drho*(Nrho-1),Nrho);
a=[0.3042,0.31231,0.342265,0.383723,0.407627,0.439499,0.336640,1.171804];
b=[0,-0.162200,-0.2372,-0.1992,-0.1992,-0.1992,0.32440,0.82755];
c=[0.811,0.250,-0.1,0,0,0,0.1735,0.2470];
rhop0=1; rhop=[0.9,0.75,0.56,0.44,0.28,1.2,2.65];
m=2; n=2;
for i=1:Nrho
    rho=rhol(i);
    if (rhop(1)<=rho && rho<=rhop(6))
        F(i)=a(1)+c(1)*(rho-rhop0)^2;
    elseif (rhop(2)<=rho && rho<=rhop(1))
        F(i)=a(2)+b(2)*(rho-rhop(1))+c(2)*(rho-rhop(1))^2;
    elseif (rhop(3)<=rho && rho<=rhop(2))
        F(i)=a(3)+b(3)*(rho-rhop(2))+c(3)*(rho-rhop(2))^2;
    elseif (rhop(4)<=rho && rho<=rhop(3))
        F(i)=a(4)+b(4)*(rho-rhop(3))+c(4)*(rho-rhop(3))^2;     
    elseif (rhop(5)<=rho && rho<=rhop(4))
        F(i)=a(5)+b(5)*(rho-rhop(4))+c(5)*(rho-rhop(4))^2;
    elseif (rho<=rhop(5))
        F(i)=(a(6)+b(6)*(rho-rhop(5))+c(6)*(rho-rhop(5)^5))*(2*rho/rhop(5)-(rho/rhop(5))^2);
    elseif (rhop(6)<=rho && rho<=rhop(7))  
        F(i)=a(7)+b(7)*(rho-rhop(6))+c(7)*(rho-rhop(6))^m;  
    elseif (rhop(7)<rho) 
        F(i)=a(8)+b(8)*(rho-rhop(7))+c(8)*(rho-rhop(7))^n;
    end
    fprintf (fid,'%12.8e\n',F(i));
end
subplot(1,4,1)
plot(rhol,F)
ylabel('\Phi(\rho) eV');
xlabel('\rho');
xlim([0 7])
grid on

%density function rho(r) (Nr values)
p=[30.134,2.1233];
psi=p(1).*exp(-p(2).*rl);
subplot(1,4,2)
plot(rl,psi)
ylabel('\psi(r); \rho(r)');
xlabel('r');
xlim([1.8 7.6])
grid on

%embeding function phi(r) (Nr values)
for i=1:length(rl)
    r=rl(i);
    if (2.15<=r && r<4.5)
        g=[0.72001055937877e1,-0.8015084e2,0.3189214778e3,-0.5754146568e3,0.421635476033e3];
        phi(i)=g(1)+g(2)/r+g(3)/r/r+g(4)/r^3+g(5)/r^4;
    elseif (4.5<=r && r<=7.62)
        g=[0.12841892e2,-0.309649556e3,0.27369890719e4,-0.10473860879e5,0.14598307517979e5];
        phi(i)=g(1)+g(2)/r+g(3)/r/r+g(4)/r^3+g(5)/r^4;
    else
        phi(i)=0.748224 + 0.63277*(2.15 - r)+ 0.72*(exp(2.96*(2.15 - r)) -1);
    end
end
subplot(1,4,3)
plot(rl,phi)
ylabel('\phi(r) eV');
xlabel('r');
xlim([1.8 7.6])
grid on
%phi=phi*2;

for i=1:length(rl)
    fprintf(fid,'%12.6e\n',psi(i));
end

Z=rl.*phi;
for i=1:length(rl)
    fprintf(fid,'%12.6e\n',Z(i));
end
subplot(1,4,4)
plot(rl,Z)
ylabel('Z=r*\rho');
xlabel('r');
fclose(fid);