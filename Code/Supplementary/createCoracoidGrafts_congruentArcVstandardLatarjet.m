function [] = createCoracoidGrafts_congruentArcVstandardLatarjet(scapulaMesh,planes,generatePlots)

%% This function serves to import in and create the base surface system for
%  running the FEA analysis of the humeral head against the glenoid.
%
%  Inputs:
%
%   scapulaMesh         surface mesh data for scapula
%   planes              planes data for scapula/humerus
%   generatePlots       flag whether to generate figures from the
%                       processing throughout the function (default = false)
%
%  Outputs
%
%   ...                 ...
%
%  This function, like parts of the main code, uses elements of the GIBBON
%  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
%  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).

    %% Check inputs

    %Participant directory
    if nargin < 2
        error('At least the scapula mesh and plane outputs from earlier functions are required as inputs.')
    end

    %Generate plots flag
    if nargin < 3
        generatePlots = false;
    else
        %Check whether it is a logical operator
        if ~islogical(generatePlots)
            error('generatePlots function inpur must be a logical of true or false')
        end
    end
    
    %% Extract coracoid from scapula
    
    %Generate settings for slicing
    snapTolerance = mean(patchEdgeLengths(scapulaMesh.scapulaF,scapulaMesh.scapulaV))/100;
    n = planes.GraftPlane.normal; %Normal direction to plane
    P = planes.GraftPlane.origin; %Point on plane

    %Slicing surface (note 3rd color data output is supressed)
    [coracoidFc,coracoidVc,~,logicSide,coracoidEc] = ...
        triSurfSlice(scapulaMesh.scapulaF,scapulaMesh.scapulaV,[],P,n,snapTolerance);

    %Visualise slice
    if generatePlots
        %Plot split planes
        cFigure; subplot(1,2,1); hold on;
        hp1 = gpatch(coracoidFc(~logicSide,:),coracoidVc,'bw','none',1);
        hp2 = gpatch(coracoidFc(logicSide,:),coracoidVc,'rw','none',1);
        legend([hp1 hp2],{'Surface above plane','Surface below plane'})
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        %Plot extracted surface and boundary
        subplot(1,2,2); hold on;
        gpatch(coracoidFc(logicSide,:),coracoidVc,'w','none',1);
        gpatch(coracoidFc(~logicSide,:),coracoidVc,'w','none',0.25);
        hp1=gpatch(coracoidEc,coracoidVc,'none','b',1,3);
        hp2=quiverVec(P,n,0.05,'k');
        legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
        axisGeom; axis manual; camlight headligth;
    end

    %Extract the faces we want to keep
    [coracoidKeepF,coracoidKeepV] = patchCleanUnused(coracoidFc(logicSide,:),coracoidVc);

    %Use the grouping function to split the extra part of the scapula that
    %comes through with the cut away
    [~,indF] = groupVertices(coracoidKeepF,coracoidKeepV,1);

    %Visualise the grouping to determine which section to grab
    %Plot ungrouped sections
    cFigure; hold on;
    %Plot grouped sections
    title('Grouped Sections. Close figure, unpause and enter number to keep...')
    gpatch(coracoidKeepF,coracoidKeepV,indF,'none');
    axisGeom;
    camlight headlight;
    colormap gjet; icolorbar;
    pause
    close
    
    %Enter the number for which group to keep
    coracoidGroup = input('Enter number of coracoid process group to keep: ');
    
    %Extract set 2, which is the glenoid
    logicKeep = logical(indF == coracoidGroup);
    [extractedCoracoidF,extractedCoracoidV] = patchCleanUnused(coracoidKeepF(logicKeep,:),coracoidKeepV);

    %Merge vertices
    [extractedCoracoidF,extractedCoracoidV] = mergeVertices(extractedCoracoidF,extractedCoracoidV);

    %Self triangulate the potentially jagged edge of the cut
    extractedCoracoidEb = patchBoundary(extractedCoracoidF,extractedCoracoidV); %Get boundary edges
    indBoundary = edgeListToCurve(extractedCoracoidEb); %Convert boundary edges to a curve list
    indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
    angleThreshold = pi*(120/180); %threshold for self triangulation
    [extractedCoracoidF,extractedCoracoidV,indBoundaryCut] = ...
        triSurfSelfTriangulateBoundary(extractedCoracoidF,extractedCoracoidV,indBoundary,angleThreshold,1);
    
    %Rotate the coracoid graft so that it sits on the global plane
    %Create the transforms for the graft
    xyPlane = createPlane([0,0,0],[1,0,0],[0,1,0]);
    graftCutPlane = createPlane(planes.GraftPlane.origin,planes.GraftPlane.normal);
    worldTransform = createBasisTransform3d(xyPlane,graftCutPlane);
    
    %Transform surface
    for pp = 1:length(extractedCoracoidV)
        extractedCoracoidV(pp,:) = transformPoint3d(extractedCoracoidV(pp,:),worldTransform);    
    end
    clear pp
    
    %Transform planes
    planes.GraftPlane.origin = transformPoint3d(planes.GraftPlane.origin,worldTransform);
    planes.GraftPlane.normal = transformVector3d(planes.GraftPlane.normal,worldTransform);
    planes.InferiorShavePlane.origin = transformPoint3d(planes.InferiorShavePlane.origin,worldTransform);
    planes.InferiorShavePlane.normal = transformVector3d(planes.InferiorShavePlane.normal,worldTransform);
    planes.AnteriorShavePlane.origin = transformPoint3d(planes.AnteriorShavePlane.origin,worldTransform);
    planes.AnteriorShavePlane.normal = transformVector3d(planes.AnteriorShavePlane.normal,worldTransform);
    
    %Now that the graft is aligned to the global plane, force the boundary
    %to have a Z-level of zero
    extractedCoracoidV(indBoundaryCut,3) = 0;

    %Visualise the boundary on the cut surface
    if generatePlots
        cFigure; hold on;
        gpatch(extractedCoracoidF,extractedCoracoidV,'bw','k');
        plotV(extractedCoracoidV(indBoundaryCut,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Boundary of Cut Coracoid Graft')
    end

    %Create a surface that closes the back of the glenoid
    %Currently uses the default 0.5 point spacing of the scapula
    [backF,backV] = regionTriMesh2D({extractedCoracoidV(indBoundaryCut,[1 2])},0.5,0,0);
    backV(:,3) = mean(extractedCoracoidV(indBoundaryCut,3)); %Add/set z-level converting to 3D mesh

    %Visualise new meshes
    if generatePlots
        cFigure; hold on;
        gpatch(extractedCoracoidF,extractedCoracoidV,'bw','k');
        gpatch(backF,backV,'gw','k');
        plotV(extractedCoracoidV(indBoundaryCut,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Filled Cut of Coracoid Graft');
    end

    %Join the two element sets
    [graftF,graftV,graftC] = joinElementSets({extractedCoracoidF,backF},...
        {extractedCoracoidV,backV});

    %Merge vertices
    [graftF,graftV] = mergeVertices(graftF,graftV);

    %Visualise the joined sets
    if generatePlots
        cFigure;
        gpatch(graftF,graftV,graftC,'k');
        colormap gjet; icolorbar;
        axisGeom;
        title('Joined Surfaces for Cut and Filled Graft');
    end

    %Check boundaries to make sure there are no holes. If there is throw an
    %error as the volumetric meshing won't work.
    if ~isempty(patchBoundary(graftF,graftV))
        error('Holes detected in graft. Stopping here as volumetric meshing won''t work');
    end
    
    %%
    
    %%%%% Inferior plane seems ok, but anterior plane seems to shift to the
    %%%%% posterior side of the graft - unclear whether this is in this
    %%%%% function or the original import surface plane???
    
    %%%%% It seems to occur in creating the first graft...
    
    disp('Stuff''s going wrong from this point on...')
    
    %% Create the standard Latarjet graft
    
    %Generate settings for slicing
    snapTolerance = mean(patchEdgeLengths(graftF,graftV))/100;
    n = planes.InferiorShavePlane.normal; %Normal direction to plane
    P = planes.InferiorShavePlane.origin; %Point on plane

    %Slicing surface (note 3rd color data output is supressed)
    [standardFc,standardVc,~,logicSide,standardEc] = ...
        triSurfSlice(graftF,graftV,[],P,n,snapTolerance);

    %Visualise slice
    if generatePlots
        %Plot split planes
        cFigure; subplot(1,2,1); hold on;
        hp1 = gpatch(standardFc(~logicSide,:),standardVc,'bw','none',1);
        hp2 = gpatch(standardFc(logicSide,:),standardVc,'rw','none',1);
        legend([hp1 hp2],{'Surface above plane','Surface below plane'})
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        %Plot extracted surface and boundary
        subplot(1,2,2); hold on;
        gpatch(standardFc(logicSide,:),standardVc,'w','none',1);
        gpatch(standardFc(~logicSide,:),standardVc,'w','none',0.25);
        hp1=gpatch(coracoidEc,standardVc,'none','b',1,3);
        hp2=quiverVec(P,n,0.05,'k');
        legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
        axisGeom; axis manual; camlight headligth;
    end

    %Extract the faces we want to keep
    %Note that these are the opposite of the logic side for this plane
    [standardKeepF,standardKeepV] = patchCleanUnused(standardFc(~logicSide,:),standardVc);

    %Merge vertices
    [extractedStandardF,extractedStandardV] = mergeVertices(standardKeepF,standardKeepV);

    %Self triangulate the potentially jagged edge of the cut
    extractedStandardEb = patchBoundary(extractedStandardF,extractedStandardV); %Get boundary edges
    indBoundary = edgeListToCurve(extractedStandardEb); %Convert boundary edges to a curve list
    indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
    angleThreshold = pi*(120/180); %threshold for self triangulation
    [extractedStandardF,extractedStandardV,indBoundaryStandard] = ...
        triSurfSelfTriangulateBoundary(extractedStandardF,extractedStandardV,indBoundary,angleThreshold,1);
    
    %Rotate the coracoid graft so that it sits on the global plane
    %Create the transforms for the graft
    standardPlane = createPlane(planes.InferiorShavePlane.origin,planes.InferiorShavePlane.normal);
    worldTransform = createBasisTransform3d(xyPlane,standardPlane);
    
    %Transform surface
    for pp = 1:length(extractedStandardV)
        extractedStandardV(pp,:) = transformPoint3d(extractedStandardV(pp,:),worldTransform);    
    end
    clear pp
    
    %Transform planes
    planes.GraftPlane.origin = transformPoint3d(planes.GraftPlane.origin,worldTransform);
    planes.GraftPlane.normal = transformVector3d(planes.GraftPlane.normal,worldTransform);
    planes.InferiorShavePlane.origin = transformPoint3d(planes.InferiorShavePlane.origin,worldTransform);
    planes.InferiorShavePlane.normal = transformVector3d(planes.InferiorShavePlane.normal,worldTransform);
    planes.AnteriorShavePlane.origin = transformPoint3d(planes.AnteriorShavePlane.origin,worldTransform);
    planes.AnteriorShavePlane.normal = transformVector3d(planes.AnteriorShavePlane.normal,worldTransform);
    
    %Now that the graft is aligned to the global plane, force the boundary
    %to have a Z-level of zero
    extractedStandardV(indBoundaryStandard,3) = 0;

    %Visualise the boundary on the cut surface
    if generatePlots
        cFigure; hold on;
        gpatch(extractedStandardF,extractedStandardV,'bw','k');
        plotV(extractedStandardV(indBoundaryStandard,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Boundary of Shaved Coracoid Graft')
    end

    %Create a surface that closes the back of the glenoid
    %Currently uses the default 0.5 point spacing of the scapula
    [standardFillF,standardFillV] = regionTriMesh2D({extractedStandardV(indBoundaryStandard,[1 2])},0.5,0,0);
    standardFillV(:,3) = mean(extractedStandardV(indBoundaryStandard,3)); %Add/set z-level converting to 3D mesh

    %Visualise new meshes
    if generatePlots
        cFigure; hold on;
        gpatch(extractedStandardF,extractedStandardV,'bw','k');
        gpatch(standardFillF,standardFillV,'gw','k');
        plotV(extractedStandardV(indBoundaryStandard,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Filled Shave of Standard Coracoid Graft');
    end

    %Join the two element sets
    [standardF,standardV,standardC] = joinElementSets({extractedStandardF,standardFillF},...
        {extractedStandardV,standardFillV});

    %Merge vertices
    [standardF,standardV] = mergeVertices(standardF,standardV);

    %Visualise the joined sets
    if generatePlots
        cFigure;
        gpatch(standardF,standardV,standardC,'k');
        colormap gjet; icolorbar;
        axisGeom;
        title('Joined Surfaces for Cut and Filled Graft');
    end

    %Check boundaries to make sure there are no holes. If there is throw an
    %error as the volumetric meshing won't work.
    if ~isempty(patchBoundary(standardF,standardV))
        error('Holes detected in graft. Stopping here as volumetric meshing won''t work');
    end
    
    %% Create the congruent arc Latarjet graft
    
    %%%%% THERE'S SOME WEIRD STUFF GOING ON WITH CREATING THE CONGRUENT
    %%%%% GRAFT HERE...UNSURE WHY THE CUT IS REMAINING SOLID...???
    
    %Generate settings for slicing
    snapTolerance = mean(patchEdgeLengths(graftF,graftV))/100;
    n = planes.AnteriorShavePlane.normal; %Normal direction to plane
    P = planes.AnteriorShavePlane.origin; %Point on plane

    %Slicing surface (note 3rd color data output is supressed)
    [congruentFc,congruentVc,~,logicSide,congruentEc] = ...
        triSurfSlice(graftF,graftV,[],P,n,snapTolerance);

    %Visualise slice
    if generatePlots
        %Plot split planes
        cFigure; subplot(1,2,1); hold on;
        hp1 = gpatch(congruentFc(~logicSide,:),congruentVc,'bw','none',1);
        hp2 = gpatch(congruentFc(logicSide,:),congruentVc,'rw','none',1);
        legend([hp1 hp2],{'Surface above plane','Surface below plane'})
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        %Plot extracted surface and boundary
        subplot(1,2,2); hold on;
        gpatch(congruentFc(logicSide,:),congruentVc,'w','none',1);
        gpatch(congruentFc(~logicSide,:),congruentVc,'w','none',0.25);
        hp1=gpatch(congruentEc,congruentVc,'none','b',1,3);
        hp2=quiverVec(P,n,0.05,'k');
        legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
        axisGeom; axis manual; camlight headligth;
    end

    %Extract the faces we want to keep
    [congruentKeepF,congruentKeepV] = patchCleanUnused(congruentFc(logicSide,:),congruentVc);

    %Merge vertices
    [extractedCongruentF,extractedCongruentV] = mergeVertices(congruentKeepF,congruentKeepV);

    %Self triangulate the potentially jagged edge of the cut
    extractedCongruentEb = patchBoundary(extractedCongruentF,extractedCongruentV); %Get boundary edges
    indBoundary = edgeListToCurve(extractedCongruentEb); %Convert boundary edges to a curve list
    indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
    angleThreshold = pi*(120/180); %threshold for self triangulation
    [extractedCongruentF,extractedCongruentV,indBoundaryCongruent] = ...
        triSurfSelfTriangulateBoundary(extractedCongruentF,extractedCongruentV,indBoundary,angleThreshold,1);
    
    %Rotate the coracoid graft so that it sits on the global plane
    %Create the transforms for the graft
    congruentPlane = createPlane(planes.AnteriorShavePlane.origin,planes.AnteriorShavePlane.normal);
    worldTransform = createBasisTransform3d(xyPlane,congruentPlane);
    
    %Transform surface
    for pp = 1:length(extractedCongruentV)
        extractedCongruentV(pp,:) = transformPoint3d(extractedCongruentV(pp,:),worldTransform);    
    end
    clear pp
    
    %Transform planes
    planes.GraftPlane.origin = transformPoint3d(planes.GraftPlane.origin,worldTransform);
    planes.GraftPlane.normal = transformVector3d(planes.GraftPlane.normal,worldTransform);
    planes.InferiorShavePlane.origin = transformPoint3d(planes.InferiorShavePlane.origin,worldTransform);
    planes.InferiorShavePlane.normal = transformVector3d(planes.InferiorShavePlane.normal,worldTransform);
    planes.AnteriorShavePlane.origin = transformPoint3d(planes.AnteriorShavePlane.origin,worldTransform);
    planes.AnteriorShavePlane.normal = transformVector3d(planes.AnteriorShavePlane.normal,worldTransform);
    
    %Now that the graft is aligned to the global plane, force the boundary
    %to have a Z-level of zero
    extractedCongruentV(indBoundaryCongruent,3) = 0;

    %Visualise the boundary on the cut surface
    if generatePlots
        cFigure; hold on;
        gpatch(extractedCongruentF,extractedCongruentV,'bw','k');
        plotV(extractedCongruentV(indBoundaryCongruent,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Boundary of Shaved Coracoid Graft')
    end

    %Create a surface that closes the back of the glenoid
    %Currently uses the default 0.5 point spacing of the scapula
    [congruentFillF,congruentFillV] = regionTriMesh2D({extractedCongruentV(indBoundaryCongruent,[1 2])},0.5,0,0);
    congruentFillV(:,3) = mean(extractedCongruentV(indBoundaryCongruent,3)); %Add/set z-level converting to 3D mesh

    %Visualise new meshes
    if generatePlots
        cFigure; hold on;
        gpatch(extractedCongruentF,extractedCongruentV,'bw','k');
        gpatch(congruentFillF,congruentFillV,'gw','k');
        plotV(extractedCongruentV(indBoundaryCongruent,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Filled Shave of Congruent Arc Coracoid Graft');
    end

    %Join the two element sets
    [congruentF,congruentV,congruentC] = joinElementSets({extractedCongruentF,congruentFillF},...
        {extractedCongruentV,congruentFillV});

    %Merge vertices
    [congruentF,congruentV] = mergeVertices(congruentF,congruentV);

    %Visualise the joined sets
    if generatePlots
        cFigure;
        gpatch(congruentF,congruentV,congruentC,'k');
        colormap gjet; icolorbar;
        axisGeom;
        title('Joined Surfaces for Cut and Filled Graft');
    end

    %Check boundaries to make sure there are no holes. If there is throw an
    %error as the volumetric meshing won't work.
    if ~isempty(patchBoundary(congruentF,congruentV))
        error('Holes detected in graft. Stopping here as volumetric meshing won''t work');
    end
    
    %%
    
    