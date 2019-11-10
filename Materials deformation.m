clc
clf
clear
%input the data
%readcsv explore the {X,Y} data to determine an appropriate range for the
data = csvread('22_0.csv',1,0);           %input the data in the csv
X=data(:,1);                               %make four vectors to store the data 
Y=data(:,2);
Z=data(:,3);
U=data(:,4);
V=data(:,5);
W=data(:,6);
x=data(:,14);
y=data(:,15);
 
u=data(:,16);
v=data(:,17);
%regular grid (xmin,xmax,ymin,ymax)        %make the suface area
Xmin=min(X);
Xmax=max(X);
Ymin=min(Y);
Ymax=max(Y);
 
xmin=min(x);
xmax=max(x);
ymin=min(y);
ymax=max(y);
inc=10;
%Creat x and y grid positions,stored as matrix variables 
[mx,my]=meshgrid(xmin:inc:xmax,ymin:inc:ymax);

% surf(mx,my);
%interpolate the displacement vector components independently onto the grid
Ui=griddata(x,y,U,mx,my,'cubic');        %the displacement vectors
Vi=griddata(x,y,V,mx,my,'cubic');
 
subplot(1,2,1)
quiver3(X,Y,Z,U,V,W);                  %plot the 3D picture
daspect ([1 1 1]);
view(0,90);
title('3D quiver plot');
xlabel('x-direction of the grid');
ylabel('y-direction of the grid');
colorbar

subplot(1,2,2)
quiver(X,Y,U,W);                       %plot the 2D picture
daspect ([1 1 1]);
view(0,90);
title('2D quiver plot');
xlabel('x-direction of the grid');
ylabel('y-direction of the grid');
colorbar

figure;
%calculate spatial derivatives of the displacement components 
[Uix,Uiy]=gradient(Ui,inc,inc);
[Vix,Viy]=gradient(Vi,inc,inc);
I=size(Uix,1);
J=size(Uiy,2);
%so we get 2D displacement gradient tensor D and Dted
for N = 1:I
  for  M = 1:J
        D=[Uix(N,M),Uiy(N,M);Vix(N,M),Viy(N,M)];
        D_t=D';                                 
        D(isnan(D))=0;               %delect the NaN data in the matrix 
        D_t(isnan(D_t))=0;
        %calculate the e and L and the shear stress
       e=(0.5.*(D+D_t));
       l=(0.5.*(D+D_t+D_t.*D));
        
       E = eig(0.5.*(D+D_t));           % create another matrix to show the tension and check the NaN value
       E_max(N,M)=max(E);               % then use the eig function to do the calculate
        
        L = eig(0.5.*(D+D_t+D_t.*D));
        L_max(N,M) = max(L);
 
    end
end
 
%plot the shear strain form e and L calculation
subplot(2,1,1);
surface(mx,my,E_max,'edgecolor','none');
colorbar
% h=surface(mx,my,E_max);
% set(h,'edgecolor','none');
daspect ([1 1 1]);
view(0,90);
title('max princple shear strain from the e');
xlabel('x-direction');
ylabel('y-direction');
colorbar

subplot(2,1,2);
surface(mx,my,L_max,'edgecolor','none');
daspect ([1 1 1]);
view(0,90);
title('max principe shear strain from the L');
xlabel('x-direction');
ylabel('y-direction');
 colorbar
 
%%%%%%%%addendum Lab1 about eig function and quiver plot
figure;
for n = 1:I
  for  data = 1:J
        D=[Uix(n,data),Uiy(n,data);Vix(n,data),Viy(n,data)];        %calculate the D again 
        D_t=D';
        D(isnan(D))=0;               %delect the NaN data in the matrix 
        D_t(isnan(D_t))=0;
        
        [vec,val] = eig(0.5.*(D+D_t));
        val2=zeros(2,1);                                            %find the largest value 
        val2(1)=val(1,1);
        val2(2)=val(2,2);
        [val2,index]=sort(val2,'descend');
        
        vec2=zeros(2,1);
        vec2(1,1)=vec(1,index(1));
        vec2(2,1)=vec(2,index(1));
        vec2(1,2)=vec(index(2),1);
        vec2(2,2)=vec(index(2),2);
        e_quiver_xmax(n,data)=val2(1)*vec2(1,1);                 %calculate the value
        e_quiver_ymax(n,data)=val2(1)*vec2(2,1);
        
        e_quiver_xmin(n,data)=val2(2)*vec2(1,2);
        e_quiver_ymin(n,data)=val2(2)*vec2(2,2);
  end
end
subplot(1,2,1)
quiver(mx,my,e_quiver_xmax,e_quiver_ymax);
daspect ([1 1 1]);
view(0,90);
title('principle vector max from e')
colorbar

subplot(1,2,2)
quiver(mx,my,e_quiver_xmin,e_quiver_ymin);
daspect ([1 1 1]);
view(0,90);
title('principle vector min from e')
colorbar

figure;
for nl = 1:I
  for  datal = 1:J
        D=[Uix(nl,datal),Uiy(nl,datal);Vix(nl,datal),Viy(nl,datal)];        %calculate the D again 
        D_t=D';
        D(isnan(D))=0;               %delect the NaN data in the matrix 
        D_t(isnan(D_t))=0;
        
        [vec_l,val_l] = eig(0.5.*(D+D_t));
        val3=zeros(2,1);                                            %find the largest value 
        val3(1)=val_l(1,1);
        val3(2)=val_l(2,2);
        [val3,index]=sort(val3,'descend');
        
        vec4=zeros(2,1);
        vec4(1,1)=vec_l(1,index(1));
        vec4(2,1)=vec_l(2,index(1));
        vec4(1,2)=vec_l(index(2),1);
        vec4(2,2)=vec_l(index(2),2);
        L_quiver_xmax(nl,datal)=val3(1)*vec4(1,1);                 %calculate the value
        L_quiver_ymax(nl,datal)=val3(1)*vec4(2,1);
        
        L_quiver_xmin(nl,datal)=val3(2)*vec4(1,2);
        L_quiver_ymin(nl,datal)=val3(2)*vec4(2,2);
  end
end
subplot(1,2,1)
quiver(mx,my,L_quiver_xmax,L_quiver_ymax);
daspect ([1 1 1]);
view(0,90);
title('principle vector max from L')
colorbar

subplot(1,2,2)
quiver(mx,my,L_quiver_xmin,L_quiver_ymin);
daspect ([1 1 1]);
view(0,90);
title('principle vector min from L')
colorbar

%%%%%%%%
Force = 275;
Area=0.000294;
load=Force/Area;

E_vector=E_max';
E_right=(E_vector(:))';

E_right(find(E_right==0))=[];

E_right1=E_right(1:2713);
Module1=load./E_right1;
mean(Module1(:))

E_right2=E_right(2714:5427);
Module2=load./E_right2;
mean(Module2(:))

E_right3=E_right(5428:8141);
Module3=load./E_right3;
mean(Module3(:))
% %%%
E_rightl=L_max;
E_rightl(find(E_rightl==0))=[];

E_rightl1=E_rightl(1:2713);
Modulel1=load./E_rightl1;
mean(Modulel1(:))

E_rightl2=E_rightl(2714:5427);
Modulel2=load./E_rightl2;
mean(Modulel2(:))

E_rightl3=E_rightl(5428:8141);
Modulel3=load./E_rightl3;
mean(Modulel3(:))
