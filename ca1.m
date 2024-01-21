clc;clear;
tic
%Import mesh data
load sp_grid.mat;
load ld_grid.mat;
load num_spgrid.mat;
load num_ldgrid.mat;
load grid_polygons.mat;
load grid_size_km.mat;
load sz.mat;
fukushima_coord = [141.0333, 37.4214];
dgrid = grid_size_km/111.32;

Se=zeros(sz+1);
%Fukushima location in the boundary coordinate system
y0=sz(1)/2 +1;
x0=sz(2)/2 +1;
lmax=0.5*sqrt(sz(1)^2+sz(2)^2);
v=ones(sz(1)+3,sz(2)+3);
for j = 1:sz(1)+3
    for i = 1:sz(2)+3
        l=sqrt((j-y0)^2+(i-x0)^2);
        v(j,i)=v(j,i)*(lmax-l)/lmax;%The farther you are from Fukushima, the slower the flow becomes
    end
end
%Diffusion coefficient
m=0.0084;
d=0.16;
w=v./max(v);
openfig('map.fig');%Show map
hold on

%Define the boundaries of the proliferation area
x=[grid_polygons(1).shape.Vertices(1,1),grid_polygons(end).shape.Vertices(3,1)];
y=[grid_polygons(1).shape.Vertices(1,2),grid_polygons(end).shape.Vertices(3,2)];
Ch = imagesc(x,y,Se);

%Show only the grid area (i.e. near Japan)
xlim([x(1)-0.5*dgrid,x(2)+0.5*dgrid]);
ylim([y(1)-0.5*dgrid,y(2)+0.5*dgrid]);

%Propose marine area nodes
sp_total=[];
for i=1:num_spgrid
    sp_total=[sp_total;sp_grid(i).diff_poly.Vertices];
end
sp_total=unique(sp_total,"rows");
%Propose all grid area nodes
gr_total=[];
for i=1:(num_spgrid+num_ldgrid)
    gr_total=[gr_total;grid_polygons(i).shape.Vertices];
end
gr_total=unique(gr_total,"rows");
%Remove the land part
for i=1:length(gr_total)
    if isempty(intersect ( find(gr_total(i,1)==sp_total(:,1)),find(  gr_total(i,2)==sp_total(:,2)) ) )
        temp1 = ceil(i/(sz(2)+1));
        temp2 = mod(i,sz(2)+1); if temp2==0 temp2=sz(2)+1; end
        Se(temp2,temp1)=NaN;
    end
end
set(Ch,'alphadata',~isnan(Se))
axis square %axis Set the axis range and aspect ratio; Square uses the same length of coordinate axis. Adjust the increments between data units accordingly

%Initialization of discharge points in the boundary coordinate system
Se=zeros(sz+1);
Se(x0-1:x0+1,y0-1:y0+1)=1;
Sd=zeros(sz(2)+3,sz(1)+3);%边界上留出2格不让改变,目的在计算邻居时防止超出索引报错
Sn=zeros(sz(2)+3,sz(1)+3);%边界上留出2格不让改变,目的在计算邻居时防止超出索引报错

%Sd is used simply to facilitate the calculation of neighbors, and the cell state is to see Se
t=0;%Iteration Initial Time
while(t<150)  %Simulate time iterations
    Sd(2:sz(1)+2,2:sz(2)+2)=Se;
   for i=2:sz(1)+2
       for j=2:sz(2)+2
%The first line is the pollutant that is transferred with the water flow, which increases and decreases; The %second and third lines are the self-diffusion of pollutants
            % flow direction is 6
%             Sd(i,j)= 0.707*v(i,j)*(Sd(i-1,j-1)+Sd(i+1,j-1)-Sd(i-1,j+1)-Sd(i+1,j+1))+v(i,j)*(Sd(i,j-1)-Sd(i,j+1)+Sd(i-1,j)+Sd(i+1,j))+...
%                     Sd(i,j)+(m*((Sd(i-1,j)-Sd(i,j))+(Sd(i+1,j)-Sd(i,j))+(1+w(i,j))*(Sd(i,j+1)-Sd(i,j))+(1-w(i,j))*(Sd(i,j-1)-Sd(i,j))))+...
%                     (m*d*((1+w(i,j))*(Sd(i-1,j+1)-Sd(i,j))+(1-w(i,j))*(Sd(i-1,j-1)-Sd(i,j))+(1+w(i,j))*(Sd(i+1,j+1)-Sd(i,j))+(1-w(i,j))*(Sd(i+1,j-1)-Sd(i,j))));
            % flow direction is 3,6,8,9
            Sd(i,j)= Sd(i,j)+(m*((Sd(i-1,j)-Sd(i,j))+(1+w(i,j))*(Sd(i+1,j)-Sd(i,j))+(1+w(i,j))*(Sd(i,j-1)-Sd(i,j))+(1-w(i,j))*(Sd(i,j+1)-Sd(i,j))))+...
                    (m*d*((1+w(i,j))*(Sd(i-1,j-1)-Sd(i,j))/d+(1-w(i,j))*(Sd(i-1,j+1)-Sd(i,j))+(1+w(i,j))*(Sd(i+1,j-1)-Sd(i,j))/d+(1-w(i,j))*(Sd(i+1,j+1)-Sd(i,j))));
            
%             Sd(i-1,j-1)=Sd(i-1,j-1)*(1-0.707*v(i,j));
%             Sd(i+1,j-1)=Sd(i+1,j-1)*(1-0.707*v(i,j));
%             Sd(i-1,j+1)=Sd(i-1,j+1)*(1+0.707*v(i,j));
%             Sd(i+1,j+1)=Sd(i+1,j+1)*(1+0.707*v(i,j));
%             Sd(i,j-1)=Sd(i,j-1)*(1-v(i,j));
%             Sd(i-1,j)=Sd(i-1,j)*(1-v(i,j));
%             Sd(i+1,j)=Sd(i+1,j)*(1-v(i,j));
%             Sd(i,j+1)=Sd(i,j+1)*(1+v(i,j));
                       

          % updates the overall cell state
            Se=Sd(2:sz(1)+2,2:sz(2)+2);
              set(Ch,"cdata",Se)
              %Uncommenting can demonstrate the process, but it will increase the time the program runs.
              %drawnow
       end
   end

    t=t+1;
                Sd(i-1-round(0.707*v(i-1,j+1)),j+1+round(0.707*v(i-1,j+1)))=Sd(i-1,j+1);
            Sd(i+1+round(0.707*v(i+1,j+1)),j+1+round(0.707*v(i+1,j+1)))=Sd(i+1,j+1);
            Sd(i+1+round(v(i+1,j)),j)=Sd(i+1,j);
            Sd(i,j+1+round(v(i,j+1)))=Sd(i,j+1);
end
toc

