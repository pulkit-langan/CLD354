syms  dr  r0 q v0 vl ts t0 c1 c2 z y  f df  tol h r

ts=117
t0=25
r0=0.1
q=64.44
v0=q/(pi*r0^2)
dr=30
dv=zeros(1,1000)
dv(1,1)=9450
tol=0.0000001
h=0.000001
c1=200  %c1=(2*k*sqrt(pi))/(rho*cp)
c2=200  %c2=n/(rho*cp)

dy = @(z,y)[(-y(1)*y(3))/y(2) ; y(3) ; y(3)^2/y(2) ; (-c1/(y(2)*sqrt(y(1))))*(y(4)-t0)+c2*(y(3)^2)]




for j=2:1000

for i=1:4
[z,yi]= ode45(dy , [0 0.5] , [ pi*r0^2 v0 dv(1,j-1) ts]);
end


[m,n]=size(yi)

f=yi(m,2)-dr*v0
for k=1:4
[z,yk]= ode45(dy , [0 0.5] , [ pi*r0^2 v0 (dv(1,j-1)+h) ts]);
end
[m,n]=size(yk)
df=((yk(m,2)-dr*v0)-f)/h 
dv(1,j)=dv(1,j-1)-f/df


if(f==0)
break
end

end

r=sqrt(yi(:,1)/pi)

plot(z,yi(:,4))
xlabel('Distance')
ylabel('Temperature')


plot(z,yi(:,3))
xlabel('Distance')
ylabel('Velocity Gradient')


plot(z,yi(:,2))
xlabel('Distance')
ylabel('velocity')


plot(z,r)
xlabel('Distance')
ylabel('Radius')
