
%Navigate to processed mesh directory
cd('..\..\Segmentation\Exported');

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

%Convert starting surface of glenoid to flat surface
surfaceV(:,3) = 0;

%% Create a sphere to map the outer curvature on to for bone
    
%Set the curvature to the radius value of 26.37 as specified by Walia
%et al. (2013), J Orthop Res, 31: 601-607 and Walia et al. (2015),
%Arthroscopy, 31: 2119-2127. Also set the minimum thickness t the
%centre of the cartilage (i.e. our coordinate system origin) to 1.14.
%But for bone...
curvatureRadius = 34.56;
baseBoneThickness = 1;

%Create a sphere using these parameters via geom3 methods. This sphere
%is positioned so that it's outer edge coincides with the minimum
%thickness away from the deep glenoid origin point.
sphere = [0 0 baseBoneThickness+curvatureRadius curvatureRadius];

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
boneBaseEb = patchBoundary(surfaceF,surfaceV);
boneBaseBC = edgeListToCurve(boneBaseEb);
boneBaseBC = boneBaseBC(1:end-1)'; %Start=End for closed curve so remove double entry

%Set start loft
startLoft = surfaceV(boneBaseBC,:);

%Create the second sketch from the end surface faces

%Get the boundary points of the two faces from the thickness approach
boneEndEb = patchBoundary(endF,endV);
boneEndBC = edgeListToCurve(boneEndEb);
boneEndBC = boneEndBC(1:end-1)'; %Start=End for closed curve so remove double entry

%Set start loft
endLoft = endV(boneEndBC,:);

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

%% Add cartilage section

%% Create a sphere to map the outer curvature on to
    
%Set the curvature to the radius value of 26.37 as specified by Walia
%et al. (2013), J Orthop Res, 31: 601-607 and Walia et al. (2015),
%Arthroscopy, 31: 2119-2127. Also set the minimum thickness t the
%centre of the cartilage (i.e. our coordinate system origin) to 1.14.
curvatureRadius = 26.37;
minThickness = 1.14;

%Create a sphere using these parameters via geom3 methods. This sphere
%is positioned so that it's outer edge coincides with the minimum
%thickness away from the deep glenoid origin point.
sphere = [0 0 minThickness+baseBoneThickness+curvatureRadius curvatureRadius];

%Loop through the vertices of the cartilage surface and create lines
%that run along the Z-axis vector. The lines structure contains the
%X,Y,Z starting points of the line followed by the X,Y,Z directions.
for vv = 1:length(endV)
    lines(vv,:) = [endV(vv,1),endV(vv,2),endV(vv,3),...
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
    gpatch(endF,endV,'b');
    gpatch(endF,pts(keepPts,:),'r');
    axisGeom;
    camlight('headlight');
    lighting phong; axis off;
end

%Create new faces and vertices objects for the end mesh
cartilageF = endF; cartilageV = pts(keepPts,:);

%% Loft between the surfaces to connect the meshes

%Create the first sketch from the base surface faces

%Get the boundary points of the two faces from the thickness approach
cartilageBaseEb = patchBoundary(endF,endV);
cartilageBaseBC = edgeListToCurve(cartilageBaseEb);
cartilageBaseBC = cartilageBaseBC(1:end-1)'; %Start=End for closed curve so remove double entry

%Set start loft
startLoft = endV(cartilageBaseBC,:);

%Create the second sketch from the end surface faces

%Get the boundary points of the two faces from the thickness approach
cartilageEndEb = patchBoundary(cartilageF,cartilageV);
cartilageEndBC = edgeListToCurve(cartilageEndEb);
cartilageEndBC = cartilageEndBC(1:end-1)'; %Start=End for closed curve so remove double entry

%Set start loft
endLoft = cartilageV(cartilageEndBC,:);

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
% cFigure;
% title('Cartilage Outer with Tri Elements','fontSize',25);
% gpatch(cartilageOuterFtri,cartilageOuterVtri,'bw','k',1);
% axisGeom;

%Join the tri element sets
Fc = {endF,cartilageF,cartilageOuterFtri};
Vc = {endV,cartilageV,cartilageOuterVtri};
[FT2,VT2,CT2] = joinElementSets(Fc,Vc);

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
    gpatch(FT2,VT2,CT2,'k',faceAlpha1);
    camlight('headlight');
    axisGeom(gca,fontSize);
    colormap(gjet(numel(Fc))); icolorbar;
end

%% Plot create surfaces together

cFigure; hold on;
gpatch(FT,VT,'bw');
gpatch(FT2,VT2,'rw');
axisGeom;

%%

%Set up export structure
stlStruct.solidNames={'boneBase'}; %names of parts
stlStruct.solidVertices={VT}; %Vertices
stlStruct.solidFaces={FT}; %Faces
stlStruct.solidNormals={[]}; %Face normals (optional)

%%%%% TODO: invert face normals here?

%Export STL file
export_STL_txt('manualBoneBase_GlenoidCoordinateSystem.stl',stlStruct);



%Set up export structure
stlStruct.solidNames={'cartilage'}; %names of parts
stlStruct.solidVertices={VT2}; %Vertices
stlStruct.solidFaces={FT2}; %Faces
stlStruct.solidNormals={[]}; %Face normals (optional)

%%%%% TODO: invert face normals here?

%Export STL file
export_STL_txt('manualCartilage_GlenoidCoordinateSystem.stl',stlStruct);

%% Test outputs from 3matic

%Load in STLs
baseSTL = import_STL('manualBoneBase.stl');
cartilageSTL = import_STL('manualCartilage.stl');

%Extract the faces and vertices
baseF = baseSTL.solidFaces{1};
baseV = baseSTL.solidVertices{1};
cartilageF = cartilageSTL.solidFaces{1};
cartilageV = cartilageSTL.solidVertices{1};

%Merge vertices
[baseF,baseV] = mergeVertices(baseF,baseV);
[cartilageF,cartilageV] = mergeVertices(cartilageF,cartilageV);

%Remesh
[baseF,baseV] = triRemeshLabel(baseF,baseV,1.0);
[cartilageF,cartilageV] = triRemeshLabel(cartilageF,cartilageV,1.0);

%Visualise
cFigure; hold on;
gpatch(baseF,baseV,'g');
gpatch(cartilageF,cartilageV,'b');
axisGeom;
camlight('headlight');
lighting phong; axis off;

%% Test humeral head cartilage code

%Import surfaces

%Load in STLs
headSTL = import_STL('IdealisedHumeralHead.stl');
surfaceSTL = import_STL('IdealisedHumeralHeadSurface.stl');

%Extract the faces and vertices
headF = headSTL.solidFaces{1};
headV = headSTL.solidVertices{1};
surfaceF = surfaceSTL.solidFaces{1};
surfaceV = surfaceSTL.solidVertices{1};

%Merge vertices
[headF,headV] = mergeVertices(headF,headV);
[surfaceF,surfaceV] = mergeVertices(surfaceF,surfaceV);

%Remesh
% [headF,headV] = triRemeshLabel(headF,headV,1.0);
% [surfaceF,surfaceV] = triRemeshLabel(surfaceF,surfaceV,1.0); %%%NOT CLOSED

%Manually create humeral head centre point
humeralCentre = [0.0116,-0.0234,25.7510];

%Manually create humeral radius
humeralRadius = 24.285;

%Shift so that head centre is at 0,0,0 to make things easier
headV = headV - humeralCentre;
surfaceV = surfaceV - humeralCentre;

%Visualise
cFigure; hold on;
gpatch(headF,headV,'g');
gpatch(surfaceF,surfaceV,'b');
axisGeom;
camlight('headlight');
lighting phong;

%Identify mean point of surface arc
meanSurfacePt = mean(surfaceV);

% % % %Visualise
% % % cFigure; hold on;
% % % gpatch(surfaceF,surfaceV,'b');
% % % scatter3(meanSurfacePt(1),meanSurfacePt(2),meanSurfacePt(3),'r');
% % % axisGeom;
% % % camlight('headlight');
% % % lighting phong; axis off;

%Identify vector from humeral centre (0,0,0) to mean surface point
surfaceToCentreVector = [0 - meanSurfacePt(1),...
    0 - meanSurfacePt(2),0 - meanSurfacePt(3)];

%Identify where this vector intersects the sphere
surfacePt = intersectLineSphere([0 0 0 surfaceToCentreVector],...
    [0 0 0 humeralRadius]);

%Take point with negative Z
intPt = surfacePt(find(surfacePt(:,3) < 0),:);

%%%% TODO: intPt doesn't quite look right???

%Visualise geom3d version of humeral sphere with above vector line
figure; hold on;
drawSphere([0 0 0 humeralRadius]);
drawLine3d([0 0 0 surfaceToCentreVector]);
drawPoint3d(intPt,'rx');
view(3);

%Visualise patch version
cFigure; hold on;
gpatch(headF,headV,'g');
gpatch(surfaceF,surfaceV,'b');
scatter3(intPt(1),intPt(2),intPt(3),'r');
axisGeom;

%Identify vector from humeral centre (0,0,0) to intersection point
surfaceToIntVector = [0 - intPt(1),0 - intPt(2),0 - intPt(3)];

%Identify surface to intersection distance
%Should just be radius actually...
surfaceToIntDist = sqrt((intPt(1) - 0)^2 + (intPt(2) - 0)^2 + (intPt(3) - 0)^2);

%Create a larger sphere using geom3d functions at the humeral centre. This
%will have the radius curvature of the Walia et al. studies of 26.85.
%This is actually probably too large as they go from 26.10 to 26.85 - we
%could probably increase by this same relative amount.
% humeralCartilageCurvatureRadius = (26.85 / 26.10) * humeralRadius;
humeralCartilageCurvatureRadius = 26.85;
bigSphere = [humeralCentre humeralCartilageCurvatureRadius];

%Create a humeral head sphere structure
humeralSphere = [humeralCentre humeralRadius];

%Calculate distance we want the big sphere to translate along the
%directional vector so that the thickness of the cartiage at the mean
%surface centre matches the 2.03 in Walia et al.
translateDist = (humeralCartilageCurvatureRadius - humeralRadius) - 2.03;

%Create proportion to multiply translation vector by
translateProportion = translateDist / surfaceToIntDist;

%Translate the sphere along the vector identified earlier
%See how much this translates without any scaling
%Plot original
figure; hold on;
drawSphere(humeralSphere,'facecolor','b');
%Translate
%Create the translation matrix
sphereTransMat = createTranslation3d((surfaceToIntVector*-1)*translateProportion);
%Shift the sphere
newCentre = transformPoint3d(humeralCentre, sphereTransMat);
%Draw the new sphere
drawSphere([newCentre humeralCartilageCurvatureRadius],'facecolor','g');
view(3);

%Loop through the vertices of the cartilage surface and create lines
%that run along the normal of the vector. To do this, we need to use the
%face normals and average those connected to the vertex
% % % cFigure; hold on;
% % % gpatch(surfaceF,surfaceV);
% % % patchNormPlot(surfaceF,surfaceV,1);
% % % axisGeom;
surfaceFnormals = patchNormal(surfaceF,surfaceV);
for vv = 1:length(surfaceV)
    %Identify faces this vertex is connected to
    indF = find(any(surfaceF == vv,2));
    %Average the face normals from the connected faces
    conFnormals = surfaceFnormals(indF,:);
    currVnormal = mean(conFnormals);
    %Add to line structure    
    lines(vv,:) = [surfaceV(vv,1),surfaceV(vv,2),surfaceV(vv,3),...
        currVnormal(1),currVnormal(2),currVnormal(3)];
    %Cleanup
    clear indF conFnormals currVnormal
end
clear vv

%Get the line to sphere intersection points using geom3d function
pts = intersectLineSphere(lines(1,:),[newCentre humeralCartilageCurvatureRadius]);

%Visualise the sphere and intersections
figure; hold on; axis equal;
drawSphere([newCentre humeralCartilageCurvatureRadius]);
drawLine3d(lines(1,:));
drawPoint3d(pts, 'rx');
view(3);

cFigure; hold on;
gpatch(surfaceF,surfaceV);
scatter3(pts(1,1),pts(1,2),pts(1,3),'r')
scatter3(pts(2,1),pts(2,2),pts(2,3),'b')
axisGeom;

%%%%% sphere intersection approach is just wrong right now...

%%%%% might just be better to do this in 3matic to create the outer
%%%%% surface...!

% % % %Test vertex normals by doing consistent addition to mesh
% % % for vv = 1:length(surfaceV)
% % %     %Identify faces this vertex is connected to
% % %     indF = find(any(surfaceF == vv,2));
% % %     %Average the face normals from the connected faces
% % %     conFnormals = surfaceFnormals(indF,:);
% % %     currVnormal = mean(conFnormals);
% % %     %Add to vertexby 3    
% % %     surfaceVadd(vv,:) = surfaceV(vv,:) + currVnormal(1,:);
% % %     %Cleanup
% % %     clear indF conFnormals currVnormal
% % % end
% % % clear vv
% % % 
% % % figure; hold on;
% % % scatter3(surfaceV(:,1),surfaceV(:,2),surfaceV(:,3));
% % % scatter3(surfaceVadd(:,1),surfaceVadd(:,2),surfaceVadd(:,3));
% % % view(3);

%%%% Vector normals look correct...


