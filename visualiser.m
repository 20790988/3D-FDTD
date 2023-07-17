clearvars -except grid source param
clf



PEC = 1;
AIR = 2;
DIELECTRIC = 3;
SUPERCONDUCTOR = 4;

%====================PLOT CONTROL====================%
to_plot = [PEC];
plot_source = true;
%===================================================%


s = size(grid);
newgrid = zeros(s);

newgrid(grid==0) = PEC;
newgrid(grid==1) = AIR;
newgrid(grid==2) = DIELECTRIC;
newgrid(grid==3) = SUPERCONDUCTOR;

temp = zeros([s(1)+2,s(2)+2, s(3)+2]);
temp(2:(s(1)+1),2:(s(2)+1),2:(s(3)+1)) = newgrid;
newgrid = temp;

material = [PEC,AIR,DIELECTRIC,SUPERCONDUCTOR];
color = {[0.5 0.5 0.5], [0.1, 1 , 1], [0.2 1 0.2 ], [0 0 1]};
figure(1);

delta_x = param.delta{1};
delta_y = param.delta{2};
delta_z = param.delta{3};

plot_x = (0:s(1)+1)*delta_x;
plot_y = (0:s(2)+1)*delta_y;
plot_z = (0:s(3)+1)*delta_z;


for i=to_plot
    
    s = size(newgrid);
    plot_grid = zeros(s);

    plot_grid(newgrid==material(i)) = 100;

    
    k = 50;
    p1 = patch(isosurface(plot_y,plot_x,plot_z,plot_grid,k));
    p2 = patch(isocaps(plot_y,plot_x,plot_z,plot_grid,k));
    
    set(p1,'FaceColor',color{i});  
    set(p2,'FaceColor',color{i});  
    set(p1,'EdgeColor','none');  
    set(p2,'EdgeColor','none');  
    
end

if plot_source
      
     s = size(newgrid);
    plot_grid = zeros(s);

    source_x = source.coord{1};
    source_y = source.coord{2};
    source_z = source.coord{3};

    plot_grid(source_x,source_y,source_z)  = 100;

    k = 50;
    p1 = patch(isosurface(plot_y,plot_x,plot_z,plot_grid,k));
    p2 = patch(isocaps(plot_y,plot_x,plot_z,plot_grid,k));

    set(p1,'FaceColor',[1 0 0]);  
    set(p2,'FaceColor',[1 0 0]);  
    set(p1,'EdgeColor','none');  
    set(p2,'EdgeColor','none');  
end

camlight('right');
lighting gouraud;
axis equal
lightangle(0,45);

p1 = patch(NaN,NaN,NaN,NaN);
p2 = patch(NaN,NaN,NaN,NaN);
p3 = patch(NaN,NaN,NaN,NaN);
p4 = patch(NaN,NaN,NaN,NaN);
p5 = patch(NaN,NaN,NaN,NaN);

set(p1,'FaceColor',color{1});  
set(p2,'FaceColor',color{2});  
set(p3,'FaceColor',color{3});  
set(p4,'FaceColor',color{4});  
set(p5,'FaceColor',[1 0 0]);  

legend([p1,p2,p3,p4,p5],{'PEC','AIR','DIELECTRIC','SUPERCONDUCTOR','SOURCE'})
xlim([0 (s(2)+2).*delta_y]);
ylim([0 (s(1)+2).*delta_x]);
zlim([0 s(3).*delta_z]);

ylabel('x (m)');
xlabel('y (m)');
zlabel('z (m)');
