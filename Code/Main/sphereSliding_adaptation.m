%%

warning off

%Navigate to surface mesh directory
cd('..\..\FEA\SurfaceMeshes');

%Load in surface file
[scapulaSTLstruct] = import_STL('Scapula.stl');

%Access the data from the STL structs

%Humeral head

%Get surface faces and vertices
scapulaF = scapulaSTLstruct.solidFaces{1}; %Faces
scapulaV = scapulaSTLstruct.solidVertices{1}; %Vertices

%Merge vertices
[scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);

%Visualise all loaded in surfaces together
cFigure;
gpatch(scapulaF,scapulaV,'kw','none');
title('Scapula');
camlight headlight;
axisGeom;

%Slice the scapula 10mm back from the glenoid origin
%Generate settings for slicing
cutLevel = -10; %Set the cut level
snapTolerance = mean(patchEdgeLengths(scapulaF,scapulaV))/100;
n = vecnormalize([0 0 -1]); %Normal direction to plane
P = [0 0 cutLevel]; %Point on plane

%Slicing surface (note 3rd color data output is supressed)
[scapulaFc,scapulaVc,~,logicSide,scapulaEc] = triSurfSlice(scapulaF,scapulaV,[],P,n,snapTolerance);

%Plot split planes
cFigure; subplot(1,2,1); hold on;
hp1 = gpatch(scapulaFc(~logicSide,:),scapulaVc,'bw','k',1);
hp2 = gpatch(scapulaFc(logicSide,:),scapulaVc,'rw','k',1);
legend([hp1 hp2],{'Surface above plane','Surface below plane'})
axisGeom; axis manual; camlight headligth;
colormap gjet;
set(gca,'FontSize',25);
%Plot extracted surface and boundary
subplot(1,2,2); hold on;
gpatch(scapulaFc(logicSide,:),scapulaVc,'w','k',1);
gpatch(scapulaFc(~logicSide,:),scapulaVc,'w','none',0.25);
hp1=gpatch(scapulaEc,scapulaVc,'none','b',1,3);
hp2=quiverVec(P,n,50,'k');
legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
axisGeom; axis manual; camlight headligth;
set(gca,'FontSize',25);

%Extract the faces we want to keep
[scapulaKeepF,scapulaKeepV] = patchCleanUnused(scapulaFc(logicSide,:),scapulaVc);

%Use the grouping function to split the extra part of the scapula that
%comes through with the cut away
[indV,indF] = groupVertices(scapulaKeepF,scapulaKeepV,1);

%Visualise the grouping
%Plot ungrouped sections
cFigure; subplot(1,2,1); hold on;
title('Ungrouped')
gpatch(scapulaKeepF,scapulaKeepV,'kw','none');
axisGeom;
camlight headlight;
%Plot grouped sections
subplot(1,2,2); hold on;
title('Grouped')
gpatch(scapulaKeepF,scapulaKeepV,indF,'none');
axisGeom;
camlight headlight;
colormap gjet; icolorbar;

%Extract set 2, which is the glenoid
logicKeep = logical(indF == 2);
[extractedGlenoidF,extractedGlenoidV] = patchCleanUnused(scapulaKeepF(logicKeep,:),scapulaKeepV);

%Merge vertices
[extractedGlenoidF,extractedGlenoidV] = mergeVertices(extractedGlenoidF,extractedGlenoidV);

%Self triangulate the potentially jagged edge of the cut
extractedGlenoidEb = patchBoundary(extractedGlenoidF,extractedGlenoidV); %Get boundary edges
indBoundary = edgeListToCurve(extractedGlenoidEb); %Convert boundary edges to a curve list
indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
angleThreshold = pi*(120/180); %threshold for self triangulation
[extractedGlenoidF,extractedGlenoidV,indBoundaryBack] = ...
    triSurfSelfTriangulateBoundary(extractedGlenoidF,extractedGlenoidV,indBoundary,angleThreshold,1);

%Force boundary to have a Z level aligned with the cut
extractedGlenoidV(indBoundaryBack,3) = cutLevel;

%Visualise the boundary on the cut surface
cFigure; hold on;
gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
camlight('headlight');
axisGeom;

%Create a surface that closes the back of the glenoid
%uses 1.5 point spacing
[backF,backV] = regionTriMesh2D({extractedGlenoidV(indBoundaryBack,[1 2])},1.5,0,0);
backV(:,3) = mean(extractedGlenoidV(indBoundaryBack,3)); %Add/set z-level converting to 3D mesh

%Visualise new meshes
cFigure; hold on;
gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
gpatch(backF,backV,'gw','k');
plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
camlight('headlight');
axisGeom;

%Join the two element sets
[glenoidF,glenoidV,glenoidC] = joinElementSets({extractedGlenoidF,backF},{extractedGlenoidV,backV});

%Merge vertices
[glenoidF,glenoidV] = mergeVertices(glenoidF,glenoidV);

%Visualise the joined sets
cFigure;
gpatch(glenoidF,glenoidV,glenoidC,'k');
colormap gjet; icolorbar;
axisGeom;

%%

%Create tetgen input structure
glenoidInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
glenoidInputStruct.Faces = glenoidF; %Boundary faces
glenoidInputStruct.Nodes = glenoidV; %Nodes of boundary
glenoidInputStruct.regionPoints= getInnerPoint(glenoidF,glenoidV); %Interior points for regions
glenoidInputStruct.holePoints = []; %Interior points for holes
glenoidInputStruct.regionA = tetVolMeanEst(glenoidF,glenoidV); %Desired tetrahedral volume for each region

%Mesh model
[glenoidMesh] = runTetGen(glenoidInputStruct); %Run tetGen

%Get outputs of mesh structure
glenoidVolE = glenoidMesh.elements; %The elements
glenoidVolV = glenoidMesh.nodes; %The vertices or nodes
glenoidVolCE = glenoidMesh.elementMaterialID; %Element material or region id
glenoidVolFb = glenoidMesh.facesBoundary; %The boundary faces
glenoidVolCb = glenoidMesh.boundaryMarker; %The boundary markers

%Visualise glenoid mesh
meshView(glenoidMesh);
title('Tetrahedral Mesh: Glenoid','FontSize',25);

%Set as V1
V1 = glenoidVolV;

%% Creating triangulated sphere surface model

[E2,V2,~] = geoSphere(3,24); %24mm radius

%Shift sphere across by it's radius and 10mm offset
contactInitialOffset = 10;
V2(:,3) = V2(:,3) + (24 + contactInitialOffset);

center_of_mass = mean(V2,1);

cFigure; hold on;
gpatch(E2,V2,'kw','k'); 
gpatch(glenoidF,glenoidVolV,glenoidC,'k');
axisGeom;

%% Joining node sets
V = [V1;V2;]; %Combined node sets
E2 = E2+size(V1,1); %Fixed element indices

% Plotting joined geometry
cFigure; hold on;
title('Joined node sets');
gpatch(glenoidVolFb,V,glenoidC,'k'); 
gpatch(E2,V,'kw','k');
axisGeom;
camlight headlight;

%% Define boundary conditions

%Supported nodes
logicRigid = glenoidC==2;
Fr = glenoidVolFb(logicRigid,:);
bcSupportList = unique(Fr(:));

%Visualize BC's
cFigure; hold on;
title('Boundary conditions model');
gpatch(glenoidVolFb,V,'kw','none',0.3); 
gpatch(E2,V,'kw','k',1); 
plotV(V(bcSupportList,:),'k.','MarkerSize',5);
axisGeom;
camlight headlight;

%% Define contact surfaces

% The rigid master surface of the sphere
F_contact_master = E2;

% The deformable slave surface of the slab
logicContactSurf1= glenoidC == 1;
%%%%needs to be flipped for some reason?
F_contact_slave = fliplr(glenoidVolFb(logicContactSurf1,:));

% Plotting surface models
cFigure; hold on;
title('Contact sets and normal directions');
gpatch(glenoidVolFb,V,'kw','none',0.3); 
hl(1) = gpatch(F_contact_master,V,'g','k',1); 
patchNormPlot(F_contact_master,V);
hl(2) = gpatch(F_contact_slave,V,'b','k',1);
patchNormPlot(F_contact_slave,V);
legend(hl,{'Master','Slave'});
axisGeom;
camlight headlight;
drawnow;

%% Defining the FEBio input structure
% See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
% manual.

%Get a template with default settings 
[febio_spec]=febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version='2.5'; 

%Module section
febio_spec.Module.ATTR.type='solid'; 

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type='static';
stepStruct.Control.time_steps=10;
stepStruct.Control.step_size=1/10;
stepStruct.Control.time_stepper.dtmin=(1/10)/100;
stepStruct.Control.time_stepper.dtmax=1/10; 
stepStruct.Control.time_stepper.max_retries=5;
stepStruct.Control.time_stepper.opt_iter=10;
stepStruct.Control.max_refs=25;
stepStruct.Control.max_ups=0;

%Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

%Step specific control section
febio_spec.Step{1}.Control=stepStruct.Control;
febio_spec.Step{1}.ATTR.id=1;
febio_spec.Step{2}.Control=stepStruct.Control;
febio_spec.Step{2}.ATTR.id=2;
    
%Material section
c1=1e-3; %Shear-modulus-like parameter
m1=8; %Material parameter setting degree of non-linearity
k_factor=1e2; %Bulk modulus factor 
k=c1*k_factor; %Bulk modulus
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=c1;
febio_spec.Material.material{1}.m1=m1;
febio_spec.Material.material{1}.c2=c1;
febio_spec.Material.material{1}.m2=-m1;
febio_spec.Material.material{1}.k=k;

febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=center_of_mass;

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name='glenoid'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(glenoidVolE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=glenoidVolE;

febio_spec.Geometry.Elements{2}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name='Sphere'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(glenoidVolE,1)+(1:1:size(E2,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=E2;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

% -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_contact_master,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=F_contact_master;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_slave,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_slave;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

%Boundary condition section 
% -> Fix boundary conditions
febio_spec.Boundary.fix{1}.ATTR.bc='x';
febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc='y';
febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc='z';
febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc=1;
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL=-8;
%just above the 10 displacement? doesn't penetrate enough with 12 & 10,000 penalty seemingly?

% % % febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat=2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='z';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc='x';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc=2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL=25; %slide displacement?

febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc='z';
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc=2;
febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL=-5; %slide displacement?

%Contact section
contactType = 'facet-to-facet-sliding';
pointSpacings = 1.0; %I think this is spacing of glenoid
switch contactType
    case 'sticky'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sticky';
        febio_spec.Contact.contact{1}.penalty=100;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.snap_tol=0;
        febio_spec.Contact.contact{1}.max_traction=0;
        febio_spec.Contact.contact{1}.search_tolerance=0.1;
    case 'facet-to-facet-sliding'
        
        %%%%steps get real small but contact works with 0.5 gaptol
        
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='facet-to-facet sliding';
        febio_spec.Contact.contact{1}.penalty=100; %increased/ worked at 10000, doesn't seem to at 1000
        febio_spec.Contact.contact{1}.auto_penalty=0; %changed
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=1;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0.5; %changed for laugon
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
    case 'sliding_with_gaps'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding_with_gaps';
        febio_spec.Contact.contact{1}.penalty=100;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.fric_coeff=0;
        febio_spec.Contact.contact{1}.fric_penalty=0;
        febio_spec.Contact.contact{1}.ktmult=1;
        febio_spec.Contact.contact{1}.seg_up=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
    case 'sliding2'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding2';
        febio_spec.Contact.contact{1}.penalty=30;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.symmetric_stiffness=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
end

% LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1; 2 1];

febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{2}.point.VAL=[0 0; 1 0; 2 1];

%Output section 

% Defining file names
febioFebFileNamePart='sphereSliding_adaptation';
febioFebFileName=[febioFebFileNamePart,'.feb']; %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='external';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

clc
[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!

%%%%%% contact parameters seemingly need to be adjusted -- penetration
%turned off auto penalty and increased penalty to 10000, struggled to
%converge -- but the contact does work...

%%

%% Import FEBio results 

if runFlag==1 %i.e. a succesful run
    
    % Importing nodal displacements from a log file
    [time_mat, N_disp_mat,~]=importFEBio_logfile(febioLogFileName_disp); %Nodal displacements    
    time_mat=[0; time_mat(:)]; %Time

    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    [CF]=vertexToFaceMeasure(glenoidVolFb,DN_magnitude);
    
    %% 
    % Plotting the simulated results using |anim8| to visualize and animate
    % deformations 
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(glenoidVolFb,V_def,CF,'k',1); %Add graphics object to animate
    hp2=gpatch(E2,V_def,'kw','none',0.3); %Add graphics object to animate
    gpatch(glenoidVolFb,V,0.5*ones(1,3),'none',0.25); %A static graphics object
    
    axisGeom; 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
        
    % Set up animation features
    animStruct.Time=time_mat; %The time vector    
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments        
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
        [CF]=vertexToFaceMeasure(glenoidVolFb,DN_magnitude); %Current color data to use
        
        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF,V_def}; %Property values for to set in order to animate
    end        
    anim8(hf,animStruct); %Initiate animation feature    
    drawnow;

end

