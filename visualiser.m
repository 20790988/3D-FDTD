clearvars -except grid
clf
% at start =
% PEC = 0
% AIR = 1
% DIEL = 2

PEC = 100;
AIR = 0;
DIELECTRIC = 50;
SUPERCONDUCTOR = 200;

newgrid = grid;
s = size(newgrid);

color = ones(s);

newgrid(newgrid==0) = PEC;
newgrid(newgrid==1) = AIR;
newgrid(newgrid==2) = DIELECTRIC;
newgrid(newgrid==3) = SUPERCONDUCTOR;

temp = zeros([s(1)+2,s(2)+2, s(3)+2]);
temp(2:(s(1)+1),2:(s(2)+1),2:(s(3)+1)) = newgrid;
newgrid = temp;

material = [SUPERCONDUCTOR,DIELECTRIC];
color = {[0.2 1 1],[0.2 1 0.2]};
figure(1);

for i=[1,2]
    
    s = size(newgrid);
    plot_grid = zeros(s);

    plot_grid(newgrid==material(i)) = 100;

    
    k = 50;
    p1 = patch(isosurface(plot_grid,k));
    p2 = patch(isocaps(plot_grid,k));
    
    set(p1,'FaceColor',color{i});  
    set(p2,'FaceColor',color{i});  
    set(p1,'EdgeColor','none');  
%     set(p2,'EdgeColor','none');  
    
end
camlight('right');
lighting gouraud;
axis equal
lightangle(0,45);

p1 = patch(NaN,NaN,NaN,NaN);
p2 = patch(NaN,NaN,NaN,NaN);
set(p1,'FaceColor',color{1});  
set(p2,'FaceColor',color{2});  

legend([p1,p2],{'Superconductor','Dielectric'})
xlim([0 s(2)+2]);
ylim([0 s(1)+2]);
zlim([0 s(3)]);
xlabel('Cells \Delta x = 50 nm')
