function createGlenoidCartilage

% @author: Aaron Fox
% Centre for Sport Research, Deakin University
% aaron.f@deakin.edu.au
% 
% This code uses the extracted surface curve of the glenoid face from
% segmented data to create a glenoid surface mesh. The surface meshes of
% the humerus and scapula in the basic glenoid coordinate system are also
% required to calculate the required thickness of the glenoid cartilage.
% The method here follows the approach specified by Buchler et al. (2002),
% whereby the minimum distance between the bony surfaces of the humeral head
% and scapula is calculated and then halved to specify a maximum constant
% thickness cartilage for the humeral head and glenoid.
%
% This function uses various aspects of the GIBBON toolbox, which obviously
% then needs to be installed for this to work.
%
% This function should be run from it's own location (i.e. Code > Main
% directory). Without this, the manoeuvering between directories will not
% be appropriate.

% References
%
% Buchler et al. (2002). A finite element model of the shoulder: Application
% to the comparison of normal and osteoarthritic joints. Clin Biomech,
% 17(9-10): 630-639.

%%%%% TODO: add description and references (i.e. Walia et al. studies) for
%%%%% the curvature method for creating glenoid cartilage

%%%%% TODO: convert this to a supplementary function as part of the main
%%%%% GIBBON FE analysis code -- include a 'method' inputting for constant
%%%%% thickness vs. curvature arc

%%%%% TODO: add reference to geom3d need...

%% Set-up

%Set starting directory
mainDir = pwd;

%Set default method
%%%%% TODO: if this isn't input as part of function
method = 'curvature'; %%% method = 'constant-thickness';

%% Load in relevant curves and meshes

%Navigate to processed mesh directory
cd('..\..\Segmentation\Exported');

%Load in STL files
%Only need the scapula and humerus if method is constant-thickness
if strcmp(method,'constant-thickness')
    
    %Import STL
    scapulaSTL = import_STL('Scapula_r_GlenoidCoordinateSystem.stl');
    humerusSTL = import_STL('Humerus_r_GlenoidCoordinateSystem.stl');
    
    %Extract the faces and vertices
    scapulaF = scapulaSTL.solidFaces{1};
    scapulaV = scapulaSTL.solidVertices{1};
    humerusF = humerusSTL.solidFaces{1};
    humerusV = humerusSTL.solidVertices{1};
    
    %Merge vertices
    [scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);
    [humerusF,humerusV] = mergeVertices(humerusF,humerusV);

end

%Load in surface
surfaceSTL = import_STL('GlenoidFace_GlenoidCoordinateSystem.stl');

%Extract the faces and vertices
surfaceF = surfaceSTL.solidFaces{1};
surfaceV = surfaceSTL.solidVertices{1};

%Merge vertices
[surfaceF,surfaceV] = mergeVertices(surfaceF,surfaceV);

%%%%% TODO: consider offsetting the cartilage surface vertices from the
%%%%% scapula so they don't sit exactly on it -- this could impact the
%%%%% contact interface

%%%%% TODO: invert face normals of the base surface so that they point back
%%%%% towards the scapula -- can do this now or later...

%% Constant thickness method
if strcmp(method,'constant-thickness')
    
    %% Identify minimum distance for cartilage thickness

    %Set a tremendously large starting minimum distance that will easily be
    %minimised by each point
    minDist = Inf;

    %Loop through the surface points
    for ss = 1:length(surfaceV)
        %Get current surface point
        sPt = surfaceV(ss,:);
        %Calculate distance to every humerus point
        currDist = sqrt((humerusV(:,1)-sPt(1)).^2+(humerusV(:,2)-sPt(2)).^2+(humerusV(:,3)-sPt(3)).^2);
        %Get minimum distance
        currMinDist = min(currDist);    
        %Replace min distance if its less than this
        if currMinDist < minDist
            minDist = currMinDist;
            %Set the current min distance points (for visualisation)
            humerusMinPt = humerusV(find(currDist == min(currMinDist),1),:);
            glenoidMinPt = sPt;
        end
        %Cleanup
        clear sPt currDist
    end
    clear ss

    %Visualise the results
    cFigure;
    title(['Minimum Distance between Surfaces: ',num2str(round(minDist,4))],'fontSize',25);
    gpatch(surfaceF,surfaceV,'b','k',0.2,1e-5);
    hold on
    gpatch(humerusF,humerusV,'b','k',0.2,1e-5);
    scatter3(humerusMinPt(:,1),humerusMinPt(:,2),humerusMinPt(:,3),50,'green','filled')
    scatter3(glenoidMinPt(:,1),glenoidMinPt(:,2),glenoidMinPt(:,3),50,'red','filled')
    axisGeom;
    axis off;

    %Set maximum cartilage thickness as half of this distance
    maxCartilageThickness = minDist/2;

    %% Create the cartilage mesh

    %Use the patchThick function to turn the surface triangles into wedges
    [cartilageE,cartilageVE,cartilageFq1,cartilageFq2] = patchThick(surfaceF,surfaceV,1,maxCartilageThickness,2);

    %Use element2patch to get patch data
    cartilageFE = element2patch(cartilageE,[],'penta6');

    %Visualise original surfaces with cartilage attached
    cFigure; hold on;
    title('Cartilage Surface on Glenoid','fontSize',25);
    gpatch(scapulaF,scapulaV,'r','k',0.2,1e-5);
    gpatch(humerusF,humerusV,'r','k',0.2,1e-5);
    gpatch(cartilageFE,cartilageVE,'bw','k',1);
    axisGeom;

    %% Create a more user-friendly mesh for the cartilage

    %Use the loft surface to create an outer layer of the glenoid cartilage

    %Create the first sketch from the base surface faces

    %Get the boundary points of the two faces from the thickness approach
    cartilageBaseEb = patchBoundary(cartilageFq1,cartilageVE);
    cartilageBaseBC = edgeListToCurve(cartilageBaseEb);
    cartilageBaseBC = cartilageBaseBC(1:end-1)'; %Start=End for closed curve so remove double entry

    %Set start loft
    startLoft = cartilageVE(cartilageBaseBC,:);

    %Create the second sketch from the end surface faces

    %Get the boundary points of the two faces from the thickness approach
    cartilageEndEb = patchBoundary(cartilageFq2,cartilageVE);
    cartilageEndBC = edgeListToCurve(cartilageEndEb);
    cartilageEndBC = cartilageEndBC(1:end-1)'; %Start=End for closed curve so remove double entry

    %Set start loft
    endLoft = cartilageVE(cartilageEndBC,:);

    %Create a guide curve (some default parameters)
    %Should just be a straight line really...
    numStepsCurve = 3; %Number of steps for the curve (creates two element height)
    p1 = mean(startLoft,1); %First point
    n1 = [0,0,1]; %First direction vector; z direction
    p2 = mean(endLoft,1); %End point
    n2 = [0,0,1]; %End direction vector; z direction
    csapsSmoothPar = 1; %Cubic smoothening spline smoothening parameter
    f = 0.05; %Extent of tangential nature to boundary curves, surface will remain approximately orthogonal to input curves for f*distance between curves
    Vg = sweepCurveSmooth(p1,p2,n1,n2,numStepsCurve,csapsSmoothPar,f);

    %Create the loft feature
    [cartilageOuterF,cartilageOuterV] = sweepLoft(startLoft,endLoft,n1,n2,Vg);

    %Visualise the loft feature
    cFigure;
    subplot(1,2,1);
    hold on;
    h(1)=plotV(Vg,'k.-','LineWidth',2);
    h(2)=plotV(startLoft,'r.-','LineWidth',2);
    h(3)=plotV(endLoft,'g.-','LineWidth',2);
    h(4)=quiverVec(p1,n1,3,'r');
    h(5)=quiverVec(p2,n2,3,'g');
    legend(h,{'Guide curve','Start section','End section','Start direction vector','End direction vector'});
    axisGeom;

    subplot(1,2,2); hold on;
    h = gpatch(cartilageOuterF,cartilageOuterV,'r','k',1);
    axisGeom;
    legend(h,{'Loften surface'});
    camlight headlight

    %Convert the outer surface from quad to tri elements to coincide with the faces elements
    [cartilageOuterFtri,cartilageOuterVtri] = quad2tri(cartilageOuterF,cartilageOuterV,'x');
    % cFigure;
    % title('Cartilage Outer with Tri Elements','fontSize',25);
    % gpatch(cartilageOuterFtri,cartilageOuterVtri,'bw','k',1);
    % axisGeom;

    %Join the tri element sets
    Fc = {cartilageFq1,cartilageFq2,cartilageOuterFtri};
    Vc = {cartilageVE,cartilageVE,cartilageOuterVtri};
    [FT,VT,CT] = joinElementSets(Fc,Vc);

    %Visualise
    fontSize = 25;
    faceAlpha1 = 1;
    cFigure;
    p=[1 3 5];
    for q=1:1:numel(Fc)
        subplot(3,2,p(q)); hold on;
        title(['Set ',num2str(q)],'FontSize',fontSize);
        gpatch(Fc{q},Vc{q},q*ones(size(Fc{q},1),1),'k',faceAlpha1);
        camlight('headlight');
        axisGeom(gca,fontSize);
        colormap(gjet(numel(Fc)));
        caxis([0.5 numel(Fc)+0.5]);
    end
    subplot(3,2,[2 4 6]); hold on;
    title('Joined sets','FontSize',fontSize);
    gpatch(FT,VT,CT,'k',faceAlpha1);
    camlight('headlight');
    axisGeom(gca,fontSize);
    colormap(gjet(numel(Fc))); icolorbar;

    %% Export glenoid cartilage as a surface model
    
    %% TODO: set as function output instead of export
    
    %% TODO: clean up mesh to be used following output (i.e. fill holes, mesh as quad etc.)
    %  This should be done in method below...

% % %     %Set up export structure
% % %     stlStruct.solidNames={'cartilage'}; %names of parts
% % %     stlStruct.solidVertices={VT}; %Vertices
% % %     stlStruct.solidFaces={FT}; %Faces
% % %     stlStruct.solidNormals={[]}; %Face normals (optional)
% % % 
% % %     %Export STL file
% % %     export_STL_txt('GlenoidCartilage_GlenoidCoordinateSystem.stl',stlStruct);
    
%% Constant thickness method
if strcmp(method,'curvature')
    
    %% Create a sphere to map the outer curvature on to
    
    %Set the curvature to the radius value of 26.37 as specified by Walia
    %et al. (2013), J Orthop Res, 31: 601-607 and Walia et al. (2015),
    %Arthroscopy, 31: 2119-2127. Also set the minimum thickness t the
    %centre of the cartilage (i.e. our coordinate system origin) to 1.14.
    curvatureRadius = 26.37;
% % %     curvatureRadius = (26.37/34.56) * 52.5457;
    minThickness = 1.14;
    
    %Create a sphere using these parameters via geom3 methods. This sphere
    %is positioned so that it's outer edge coincides with the minimum
    %thickness away from the deep glenoid origin point.
    sphere = [0 0 minThickness+curvatureRadius curvatureRadius];
    
    %Loop through the vertices of the cartilage surface and create lines
    %that run along the Z-axis vector. The lines structure contains the
    %X,Y,Z starting points of the line followed by the X,Y,Z directions.
    for vv = 1:length(surfaceV)
        lines(vv,:) = [surfaceV(vv,1),surfaceV(vv,2),surfaceV(vv,3),...
            0,0,1];
    end
    clear vv
    
    %Get the line to sphere intersection points using geom3d function
    pts = intersectLineSphere(lines,sphere);
    
    %Find the first points by identifying those that the z-coordinates are to
    %the right (i.e. smaller number) than the centre - we can use the radius
    %given how close it is to the z-coordinate origin
    keepPts = find(pts(:,3) < 26.37);
    
    %TODO: set a 'generatePlots' boolean function input
    if generatePlots
        %Visualise the sphere and intersections
        figure; hold on; axis equal;
        drawSphere(sphere);
        drawLine3d(lines);
        drawPoint3d(pts(keepPts,:), 'rx');
        view(3);

        %Plot the two surfaces as patch data
        cFigure; hold on
        title('Glenoid Surface & Sphere Projection','fontSize',20);
        gpatch(surfaceF,surfaceV,'b');
        gpatch(surfaceF,pts(keepPts,:),'r');
        axisGeom;
        camlight('headlight');
        lighting phong; axis off;
    end
    
    %Create new faces and vertices objects for the end mesh
    endF = surfaceF; endV = pts(keepPts,:);
    
    %% Loft between the surfaces to connect the meshes

    %Create the first sketch from the base surface faces

    %Get the boundary points of the two faces from the thickness approach
    cartilageBaseEb = patchBoundary(surfaceF,surfaceV);
    cartilageBaseBC = edgeListToCurve(cartilageBaseEb);
    cartilageBaseBC = cartilageBaseBC(1:end-1)'; %Start=End for closed curve so remove double entry

    %Set start loft
    startLoft = surfaceV(cartilageBaseBC,:);

    %Create the second sketch from the end surface faces

    %Get the boundary points of the two faces from the thickness approach
    cartilageEndEb = patchBoundary(endF,endV);
    cartilageEndBC = edgeListToCurve(cartilageEndEb);
    cartilageEndBC = cartilageEndBC(1:end-1)'; %Start=End for closed curve so remove double entry

    %Set start loft
    endLoft = endV(cartilageEndBC,:);

    %Create a guide curve (some default parameters)
    %Should just be a straight line really...
    numStepsCurve = 3; %Number of steps for the curve (creates two element height)
    p1 = mean(startLoft,1); %First point
    n1 = [0,0,1]; %First direction vector; z direction
    p2 = mean(endLoft,1); %End point
    n2 = [0,0,1]; %End direction vector; z direction
    csapsSmoothPar = 1; %Cubic smoothening spline smoothening parameter
    f = 0.05; %Extent of tangential nature to boundary curves, surface will remain approximately orthogonal to input curves for f*distance between curves
    Vg = sweepCurveSmooth(p1,p2,n1,n2,numStepsCurve,csapsSmoothPar,f);

    %Create the loft feature
    [outerF,outerV] = sweepLoft(startLoft,endLoft,n1,n2,Vg);

    %Visualise the loft feature
    if generatePlots
        cFigure;
        subplot(1,2,1);
        hold on;
        h(1)=plotV(Vg,'k.-','LineWidth',2);
        h(2)=plotV(startLoft,'r.-','LineWidth',2);
        h(3)=plotV(endLoft,'g.-','LineWidth',2);
        h(4)=quiverVec(p1,n1,3,'r');
        h(5)=quiverVec(p2,n2,3,'g');
        legend(h,{'Guide curve','Start section','End section','Start direction vector','End direction vector'});
        axisGeom;

        subplot(1,2,2); hold on;
        h = gpatch(outerF,outerV,'r','k',1);
        axisGeom;
        legend(h,{'Loften surface'});
        camlight headlight
    end

    %Convert the outer surface from quad to tri elements to coincide with the faces elements
    [outerFtri,outerVtri] = quad2tri(outerF,outerV,'x');
    % cFigure;
    % title('Cartilage Outer with Tri Elements','fontSize',25);
    % gpatch(cartilageOuterFtri,cartilageOuterVtri,'bw','k',1);
    % axisGeom;

    %Join the tri element sets
    Fc = {surfaceF,endF,outerFtri};
    Vc = {surfaceV,endV,outerVtri};
    [FT,VT,CT] = joinElementSets(Fc,Vc);

    %Visualise
    if generatePlots
        fontSize = 25;
        faceAlpha1 = 1;
        cFigure;
        p=[1 3 5];
        for q=1:1:numel(Fc)
            subplot(3,2,p(q)); hold on;
            title(['Set ',num2str(q)],'FontSize',fontSize);
            gpatch(Fc{q},Vc{q},q*ones(size(Fc{q},1),1),'k',faceAlpha1);
            camlight('headlight');
            axisGeom(gca,fontSize);
            colormap(gjet(numel(Fc)));
            caxis([0.5 numel(Fc)+0.5]);
        end
        subplot(3,2,[2 4 6]); hold on;
        title('Joined sets','FontSize',fontSize);
        gpatch(FT,VT,CT,'k',faceAlpha1);
        camlight('headlight');
        axisGeom(gca,fontSize);
        colormap(gjet(numel(Fc))); icolorbar;
    end

    %% Clean up the joined element sets
    
    %%%%% Seems difficult to cleanup in Matlab with GIBBON -- might still
    %%%%% need to export and fix up in MeshLab...
    
% % %     %Search for and cleanup any holes in generated surface mesh
% % %     [Fclosed,Vclosed] = triSurfCloseHoles(FT,VT);
% % %     
% % %     
% % %     [Fre,Vre] = triRemeshLabel(Fclosed,Vclosed,0.5);

    %% Export glenoid cartilage as a surface model
    
    %% TODO: set as function output instead of export
    
    %% TODO: clean up mesh to be used following output (i.e. fill holes, mesh as quad etc.)
    %  This should be done in method below...

    %Set up export structure
    stlStruct.solidNames={'cartilage'}; %names of parts
    stlStruct.solidVertices={VT}; %Vertices
    stlStruct.solidFaces={FT}; %Faces
    stlStruct.solidNormals={[]}; %Face normals (optional)
    
    %%%%% TODO: invert face normals here?

    %Export STL file
    export_STL_txt('GlenoidCartilage_CurvatureMethod_GlenoidCoordinateSystem.stl',stlStruct);
    
end
    
end

%% ----- End of createGlenoidCartilage.m ----- %%