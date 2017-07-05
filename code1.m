syms r a b z q vr y c1 c2 y0 h cp p dy dv

q=64.44

vr=117340

z=zeros(41,1)

for i=2:41
   z(i,1)=z(i-1,1)+0.0125
end

r=zeros(41,1)

r(1,1)=0.1
r(41,1)=sqrt(q/(vr*pi))

b=r(1,1)
a=(r(41,1)-b)/z(41,1)

clear r,z

c1=-1
c2=0.001
y0=25



dy= @ (z,y)[ -y(1)*a/(a*z+b) ; c1*(y(2)-y0)/(a*z+b) + c2*y(1)/(a*z+b)^6 ];

for i=1:2
[z,yi]=ode45(dy,[0 0.5],[0.1 100]);
end

dv=-2*a/a*z+b

plot(dv,yi(:,1),'g')
xlabel('dv/dz') 
ylabel('Viscocity') 

plot(z,a*z+b,'r')
xlabel('Distance-z')
ylabel('Dimensionless Diameter- D')

plot(z,yi(:,1),'b')
xlabel('Distance-z')
ylabel('Viscocity')

c1=-2*h/p*cp
c2=4*q^2*a^2/pi^2*p*cp