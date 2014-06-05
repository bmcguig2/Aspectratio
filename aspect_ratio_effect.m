 %include path of comsol functions
path(path,'C:\Program Files\COMSOL\COMSOL43\petsc');
path(path,'C:\Program Files\COMSOL\COMSOL43\mli')

% load the model
disp('Loading Model')
model=mphload('C:\Users\bmcgu_000\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\aspect_ratio_effect.mph');


%Change quantum dot aspect ratio
disp('Changing Geometry')
model.geom('geom1').feature('elp1').setIndex('semiaxes', '6', 0);
model.geom('geom1').feature('elp1').setIndex('semiaxes', '8', 1);
model.geom('geom1').runAll;

%Remesh 
disp('Remeshing')
model.mesh('mesh4').run; %Remesh lattice strain mesh
model.mesh('mesh5').run; %Remesh valence/conduction band mesh

%Solve lattice strain
disp('Solving Lattice Strain')
model.sol('sol7').runAll;

% Make and save Valence matrices
disp('Making Valence Stiffness Matrix')
K=mphmatrix(model,'sol4','out','K'); %load Stiffness matrix for valence band
disp('Saving Stiffness Matrix')
PetscBinaryWrite('C:\Users\bmcgu_000\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\Valence\K.petsc',K.K);
clear K

disp('Making Valence Mass Matrix')
D=mphmatrix(model,'sol4','out','D'); % load Mass Matrix for valence band
disp('Saving Stiffness Matrix')
PetscBinaryWrite('C:\Users\bmcgu_000\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\Valence\D.petsc',D.D);
clear D

% Make and save Conduction matrices
disp('Making Conduction Stiffness Matrix')
K=mphmatrix(model,'sol2','out','K'); %load Stiffness matrix for valence band
disp('Saving Stiffness Matrix')
PetscBinaryWrite('C:\Users\bmcgu_000\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\Conduction\K.petsc',K.K);
clear K

disp('Making Conduction Mass Matrix')
D=mphmatrix(model,'sol2','out','D'); % load Mass Matrix for valence band
disp('Saving Stiffness Matrix')
PetscBinaryWrite('C:\Users\bmcgu_000\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\Conduction\D.petsc',D.D);
clear D



%eigen=PetscBinaryRead('C:\Users\Brian\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\Valence\eigenvals.out');


%extract isosurface solution data from comsol and plot in matlab with
%transparency

% dat=mphplot(model,'pg30');
% mphplot(dat);
% alpha(0.4);


%retreive mesh info
 info=mphxmeshinfo(model,'soltag','sol4');
 %solCoords = info.dofs.coords;
 
 %plot nodal solution on extracted mesh nodes
  %scatter3(info.nodes.coords(1,:),info.nodes.coords(2,:),info.nodes.coords(3,:),20,V(:,1),'o','filled')
  
 % test=mpheval(V(:,1))
  
  [V,ei]=eigs(K.K,D.D,1,0.8*1.602e-19);
  
 
  mask=(info.dofs.nameinds==4);
  mapping=info.dofs.nodes(mask,1);
  mapping=mapping+1;
  x=info.nodes.coords(1,mapping)*1e9;
  y=info.nodes.coords(2,mapping)*1e9;
  z=info.nodes.coords(3,mapping)*1e9;
  scatter3(x,y,z,10,V(:,1)','o','filled')
  hold on
  mphgeom(model,'geom1','Facealpha',0.1)

  
  %scatter3(info.nodes.coords(1,mapping),info.nodes.coords(2,mapping),info.nodes.coords(3,mapping),20,V(:,1),'o','filled')
  
  %MAke my own isosurface plotter
  z=linspace(min(info.nodes.coords(3,mapping)),max(info.nodes.coords(3,mapping)),150);
  y=linspace(min(info.nodes.coords(2,mapping)),max(info.nodes.coords(2,mapping)),150);
  x=linspace(min(info.nodes.coords(1,mapping)),max(info.nodes.coords(1,mapping)),150);
  [X,Y,Z]=ndgrid(x,y,z);
  [test]=griddata(info.nodes.coords(1,mapping),info.nodes.coords(2,mapping),info.nodes.coords(3,mapping),V(:,1),X,Y,Z);
  isosurface(X,Y,Z,test,-0.09)
  
   in=PetscBinaryRead('C:\Users\bmcgu_000\Desktop\Research\Combined APT and XSTM Paper\Simulations\Aspect_ratio_effect\eigenvalue.petsc','cell',10000);

  