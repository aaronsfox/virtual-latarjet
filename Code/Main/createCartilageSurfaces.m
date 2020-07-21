function createCartilageSurfaces(method,generatePlots)

% @author: Aaron Fox
% Centre for Sport Research, Deakin University
% aaron.f@deakin.edu.au
%
% Inputs:
%   method:         string of 'constant-thickness' or 'curvature' (default)
%                   dictating the method to use in creating cartilage surfaces
%   generatePlots:  logical as to whether to visualise results as the
%                   function runs
% 
% This code uses the extracted surface curve of the glenoid face from
% segmented data to create a glenoid surface mesh. The surface meshes of
% the humerus and scapula in the basic glenoid coordinate system are also
% required to calculate the required thickness of the glenoid cartilage if
% using the constant thickness method. The function also calculates
% corresponding humeral head surface cartilage based on an idealised cut
% 'dome' shape.
%
% Constant Thickness Method:
% This method here follows the approach specified by Buchler et al. (2002),
% whereby the minimum distance between the bony surfaces of the humeral head
% and scapula is calculated and then halved to specify a maximum constant
% thickness cartilage for the humeral head and glenoid.
%
% Curvature Method:
% This method uses the approach outlined in two studies by Walia et al.
% (2013,2015) and considers the minimum thickness of the cartilage surfaces
% at the central point and then calculates the rest of the surface based on
% the radius of curvature of the object.
%
% This function uses various aspects of the GIBBON toolbox, found at:
% http://gibboncode.org/, which obviously then needs to be installed for
% this function to work.
%
% This function also uses elements of the geom3d package, found at:
% https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d, which
% again, also needs to be installed for this to work.
%
% This function should be run from it's own location (i.e. Code > Main
% directory). Without this, the manoeuvering between directories will not
% be appropriate.

% References
%
% Buchler et al. (2002). A finite element model of the shoulder: Application
% to the comparison of normal and osteoarthritic joints. Clin Biomech,
% 17: 630-639.
%
% Walia et al. (2013). Theoretical model of the effect of combined
% glenohumeral bone defects on anterior shoulder instability: A finite
% element approach. J Orthop Res, 31: 601-607.
%
% Walia et al. (2015). Influence of combined Hill-Sachs and bony Bankart
% defects on range of motion in anterior instability in a finite element
% model, Arthroscopy, 31: 2119-2127.

%% Set-up

    %Set default method if necessary, else check input
    if nargin < 1
        method = 'curvature';
    else
        %Check inputs
        if ~any(strcmp({'curvature','constant-thickness'},method))
            error('Method input must be curvature or constant-thickness')
        end
    end
    
    %Set default for generate plots and check inputs
    if nargin < 2
        generatePlots = false;
    else
        %Check inputs
        if ~islogical(generatePlots)
            error('Generate plots input must be a logical of true or false')
        end
    end
    
    %% Load in relevant curves and meshes

    %Navigate to processed mesh directory
    cd('..\..\Segmentation\Exported');

    %Load in STL files
    %Only need the scapula and humerus if method is constant-thickness
    %Only need the outer humeral surface if using curvature method
    if strcmp(method,'constant-thickness')

        %Import STL
        scapulaSTL = import_STL('Scapula_r.stl');
        humerusSTL = import_STL('Humerus_r.stl');

        %Extract the faces and vertices
        scapulaF = scapulaSTL.solidFaces{1};
        scapulaV = scapulaSTL.solidVertices{1};
        humerusF = humerusSTL.solidFaces{1};
        humerusV = humerusSTL.solidVertices{1};

        %Merge vertices
        [scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);
        [humerusF,humerusV] = mergeVertices(humerusF,humerusV);
        
    elseif strcmp(method,'curvature')
        
        %Import STL
        outerHeadSTL = import_STL('HumeralCartilageOuter.stl');

        %Extract the faces and vertices
        outerHeadF = outerHeadSTL.solidFaces{1};
        outerHeadV = outerHeadSTL.solidVertices{1};

        %Merge vertices
        [outerHeadF,outerHeadV] = mergeVertices(outerHeadF,outerHeadV);

    end

    %Load in surfaces
    surfaceSTL = import_STL('GlenoidFace.stl');
    headSTL = import_STL('HumeralCartilageBase.stl');
    
    %Extract the faces and vertices
    surfaceF = surfaceSTL.solidFaces{1};
    surfaceV = surfaceSTL.solidVertices{1};
    headF = headSTL.solidFaces{1};
    headV = headSTL.solidVertices{1};

    %Merge vertices
    [surfaceF,surfaceV] = mergeVertices(surfaceF,surfaceV);
    [headF,headV] = mergeVertices(headF,headV);
   
    %%%%% TODO: consider offsetting the cartilage surface vertices from the
    %%%%% scapula so they don't sit exactly on it -- this could impact the
    %%%%% contact interface

    %% Constant thickness method
    if strcmp(method,'constant-thickness')
        %% CREATE GLENOID CARTILAGE SURFACE MESH
        
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
        if generatePlots
            cFigure;
            title(['Minimum Distance between Surfaces: ',num2str(round(minDist,4))],'fontSize',25);
            gpatch(surfaceF,surfaceV,'b','k',0.2,1e-5);
            hold on
            gpatch(humerusF,humerusV,'b','k',0.2,1e-5);
            scatter3(humerusMinPt(:,1),humerusMinPt(:,2),humerusMinPt(:,3),50,'green','filled')
            scatter3(glenoidMinPt(:,1),glenoidMinPt(:,2),glenoidMinPt(:,3),50,'red','filled')
            axisGeom;
            axis off;
        end

        %Set maximum cartilage thickness as half of this distance
        maxCartilageThickness = minDist/2;

        %% Create the cartilage mesh

        %Use the patchThick function to turn the surface triangles into wedges
        [cartilageE,cartilageVE,cartilageFq1,cartilageFq2] = patchThick(surfaceF,surfaceV,1,maxCartilageThickness,2);

        %Use element2patch to get patch data
        cartilageFE = element2patch(cartilageE,[],'penta6');

        %Visualise original surfaces with cartilage attached
        if generatePlots
            cFigure; hold on;
            title('Cartilage Surface on Glenoid','fontSize',25);
            gpatch(scapulaF,scapulaV,'r','k',0.2,1e-5);
            gpatch(humerusF,humerusV,'r','k',0.2,1e-5);
            gpatch(cartilageFE,cartilageVE,'bw','k',1);
            axisGeom;
        end

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
            h = gpatch(cartilageOuterF,cartilageOuterV,'r','k',1);
            axisGeom;
            legend(h,{'Loften surface'});
            camlight headlight
        end

        %Convert the outer surface from quad to tri elements to coincide with the faces elements
        [cartilageOuterFtri,cartilageOuterVtri] = quad2tri(cartilageOuterF,cartilageOuterV,'x');

        %Join the tri element sets
        Fc = {cartilageFq1,cartilageFq2,cartilageOuterFtri};
        Vc = {cartilageVE,cartilageVE,cartilageOuterVtri};
        [FT,VT,CT] = joinElementSets(Fc,Vc);

        %Visualise
        if generatePlots
            cFigure;
            p=[1 3 5];
            for q=1:1:numel(Fc)
                subplot(3,2,p(q)); hold on;
                title(['Set ',num2str(q)],'FontSize',25);
                gpatch(Fc{q},Vc{q},q*ones(size(Fc{q},1),1),'k',1);
                camlight('headlight');
                axisGeom(gca,25);
                colormap(gjet(numel(Fc)));
                caxis([0.5 numel(Fc)+0.5]);
            end
            subplot(3,2,[2 4 6]); hold on;
            title('Joined sets','FontSize',25);
            gpatch(FT,VT,CT,'k',1);
            camlight('headlight');
            axisGeom(gca,25);
            colormap(gjet(numel(Fc))); icolorbar;
        end

        %% Export glenoid cartilage as a surface model

        %Set up export structure
        stlStruct.solidNames={'glenoidCartilage'}; %names of parts
        stlStruct.solidVertices={VT}; %Vertices
        stlStruct.solidFaces={FT}; %Faces
        stlStruct.solidNormals={[]}; %Face normals (optional)

        %Export STL file
        export_STL_txt('..\Processed\GlenoidCartilage_ConstantThickness.stl',stlStruct);
        
        %% CREATE HUMERAL CARTILAGE
        
        %% Create the cartilage mesh

        %Use the patchThick function to turn the surface triangles into wedges
        [cartilageE,cartilageVE,cartilageFq1,cartilageFq2] = patchThick(headF,headV,1,maxCartilageThickness,2);

        %Use element2patch to get patch data
        cartilageFE = element2patch(cartilageE,[],'penta6');

        %Visualise original surfaces with cartilage attached
        %Note that this won't match up on the humerus because this
        %cartilage surface is from an idealised sphere
        if generatePlots
            cFigure; hold on;
            title('Cartilage Surface on Humerus','fontSize',25);
            gpatch(scapulaF,scapulaV,'r','k',0.2,1e-5);
            gpatch(humerusF,humerusV,'r','k',0.2,1e-5);
            gpatch(cartilageFE,cartilageVE,'bw','k',1);
            axisGeom;
        end
        
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
            h = gpatch(cartilageOuterF,cartilageOuterV,'r','k',1);
            axisGeom;
            legend(h,{'Loften surface'});
            camlight headlight
        end

        %Convert the outer surface from quad to tri elements to coincide with the faces elements
        [cartilageOuterFtri,cartilageOuterVtri] = quad2tri(cartilageOuterF,cartilageOuterV,'x');

        %Join the tri element sets
        Fc = {cartilageFq1,cartilageFq2,cartilageOuterFtri};
        Vc = {cartilageVE,cartilageVE,cartilageOuterVtri};
        [FT,VT,CT] = joinElementSets(Fc,Vc);

        %Visualise
        if generatePlots
            cFigure;
            p=[1 3 5];
            for q=1:1:numel(Fc)
                subplot(3,2,p(q)); hold on;
                title(['Set ',num2str(q)],'FontSize',25);
                gpatch(Fc{q},Vc{q},q*ones(size(Fc{q},1),1),'k',1);
                camlight('headlight');
                axisGeom(gca,25);
                colormap(gjet(numel(Fc)));
                caxis([0.5 numel(Fc)+0.5]);
            end
            subplot(3,2,[2 4 6]); hold on;
            title('Joined sets','FontSize',25);
            gpatch(FT,VT,CT,'k',1);
            camlight('headlight');
            axisGeom(gca,25);
            colormap(gjet(numel(Fc))); icolorbar;
        end

        %% Export humeral head cartilage as a surface model

        %Set up export structure
        stlStruct.solidNames={'humeralCartilage'}; %names of parts
        stlStruct.solidVertices={VT}; %Vertices
        stlStruct.solidFaces={FT}; %Faces
        stlStruct.solidNormals={[]}; %Face normals (optional)

        %Export STL file
        export_STL_txt('..\Processed\HumeralCartilage_ConstantThickness.stl',stlStruct);
    
    %% Constant thickness method
    elseif strcmp(method,'curvature')
        
        %% CREATE GLENOID CARTILAGE MESH
        
        %% Create a sphere to map the outer curvature on to

        %Set the curvature of the outer surface of the cartilage to a value
        %relative to the inner 'bone' surface of the glenoid. Walia et al.
        %(2013), J Orthop Res, 31: 601-607 and Walia et al. (2015),
        %Arthroscopy, 31: 2119-2127 use values of 26.37 and 34.56 for
        %curvature of the glenoid cartilage and bone, respectively. We can
        %use the ratio of these, relative to our curvature radius of the
        %glenoid bone surface of 52.5457 to create a similarly structured
        %curvature between bone and cartilage. We also set the minimum
        %thickness at the centre of the cartilage to be 1.14 as per these
        %studies.
        curvatureRadius = (26.37/34.56) * 52.5457;
        minThickness = 1.14;

        %Create a sphere using these parameters via geom3d methods. This sphere
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
        keepPts = find(pts(:,3) < curvatureRadius);

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

        %Join the tri element sets
        Fc = {surfaceF,endF,outerFtri};
        Vc = {surfaceV,endV,outerVtri};
        [FT,VT,CT] = joinElementSets(Fc,Vc);

        %Visualise
        if generatePlots
            cFigure;
            p=[1 3 5];
            for q=1:1:numel(Fc)
                subplot(3,2,p(q)); hold on;
                title(['Set ',num2str(q)],'FontSize',25);
                gpatch(Fc{q},Vc{q},q*ones(size(Fc{q},1),1),'k',1);
                camlight('headlight');
                axisGeom(gca,25);
                colormap(gjet(numel(Fc)));
                caxis([0.5 numel(Fc)+0.5]);
            end
            subplot(3,2,[2 4 6]); hold on;
            title('Joined sets','FontSize',25);
            gpatch(FT,VT,CT,'k',1);
            camlight('headlight');
            axisGeom(gca,25);
            colormap(gjet(numel(Fc))); icolorbar;
        end

        %% Export glenoid cartilage as a surface model

        %Set up export structure
        stlStruct.solidNames={'glenoidCartilage'}; %names of parts
        stlStruct.solidVertices={VT}; %Vertices
        stlStruct.solidFaces={FT}; %Faces
        stlStruct.solidNormals={[]}; %Face normals (optional)

        %Export STL file
        export_STL_txt('..\Processed\GlenoidCartilage_Curvature.stl',stlStruct);
        
        %% CREATE HUMERAL HEAD CARTILAGE MESH
        
        %% Loft between the humeral head surfaces to connect the meshes

        %Create the first sketch from the base surface faces

        %Get the boundary points of the two faces from the thickness approach
        cartilageBaseEb = patchBoundary(headF,headV);
        cartilageBaseBC = edgeListToCurve(cartilageBaseEb);
        cartilageBaseBC = cartilageBaseBC(1:end-1)'; %Start=End for closed curve so remove double entry

        %Set start loft
        startLoft = headV(cartilageBaseBC,:);

        %Create the second sketch from the end surface faces

        %Get the boundary points of the two faces from the thickness approach
        cartilageEndEb = patchBoundary(outerHeadF,outerHeadV);
        cartilageEndBC = edgeListToCurve(cartilageEndEb);
        cartilageEndBC = cartilageEndBC(1:end-1)'; %Start=End for closed curve so remove double entry

        %Set start loft
        endLoft = outerHeadV(cartilageEndBC,:);
        
        %Interpolate end loft to same number of points as start loft
        x = 1:1:length(endLoft);
        xq = 1:length(endLoft)/length(startLoft):length(endLoft);
        endLoft = interp1(x,endLoft,xq);

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

        %Join the tri element sets
        Fc = {headF,outerHeadF,outerFtri};
        Vc = {headV,outerHeadV,outerVtri};
        [FT,VT,CT] = joinElementSets(Fc,Vc);

        %Visualise
        if generatePlots
            cFigure;
            p=[1 3 5];
            for q=1:1:numel(Fc)
                subplot(3,2,p(q)); hold on;
                title(['Set ',num2str(q)],'FontSize',25);
                gpatch(Fc{q},Vc{q},q*ones(size(Fc{q},1),1),'k',1);
                camlight('headlight');
                axisGeom(gca,25);
                colormap(gjet(numel(Fc)));
                caxis([0.5 numel(Fc)+0.5]);
            end
            subplot(3,2,[2 4 6]); hold on;
            title('Joined sets','FontSize',25);
            gpatch(FT,VT,CT,'k',1);
            camlight('headlight');
            axisGeom(gca,25);
            colormap(gjet(numel(Fc))); icolorbar;
        end

        %% Export glenoid cartilage as a surface model

        %Set up export structure
        stlStruct.solidNames={'humeralCartilage'}; %names of parts
        stlStruct.solidVertices={VT}; %Vertices
        stlStruct.solidFaces={FT}; %Faces
        stlStruct.solidNormals={[]}; %Face normals (optional)

        %Export STL file
        export_STL_txt('..\Processed\HumeralCartilage_Curvature.stl',stlStruct);

    end
    
end

%% ----- End of createCartilageSurfaces.m ----- %%