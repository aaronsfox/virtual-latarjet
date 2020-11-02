function [glenoidDefectMesh] = createAnteriorBankartDefects(glenoidMesh,landmarks,pointSpacing,generatePlots)
    
    %% TODO: consider inlcluding Latarjet process with this? Or separate

    %% This function serves to import in and create the base surface system for
    %  running the FEA analysis of the humeral head against the glenoid.
    %
    %  Inputs:
    %
    %   glenoidMesh         structure of glenoid mesh containing faces and
    %                       vertices as 'glenoidF' and 'glenoidV'
    %   landmarks           structure containing glenoid landmarks 
    %   pointSpacing        desired point spacing for remeshing if
    %                       required. Defaults to empty if no input which
    %                       means no remeshing is performed.
    %   generatePlots       flag whether to generate figures from the
    %                       processing throughout the function (default = false)
    %
    %  Outputs
    %
    %   glenoidDefectMesh   structure containing the glenoid surface meshes
    %                       (faces and vertices) of the different defect
    %                       sizes. These no longer have the coracoid aspect
    %                       as it is removed during the defect creation.
    %
    %  This function, like parts of the main code, uses elements of the GIBBON
    %  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
    %  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).
    %
    %  TODO: also add reference to geom2d
    %  (https://www.mathworks.com/matlabcentral/fileexchange/7844-geom2d)
    %  
    %  References
    %
    %  Klemt et al. (2019). The critical size of a defect in the glenoid
    %  causing anterior instability of the shoulder after a Bankart repair,
    %  under physiological joint loading. Bone Joint J, 101-B: 68-74.
    
    %% Check inputs

    %Participant directory
    if nargin < 2
        error('Structures containing the glenoid mesh as well as landmarks are required')
    end
    
    %Remeshing
    if nargin < 3
        pointSpacing = [];
    end
    
    %Generate plots flag
    if nargin < 4
        generatePlots = false;
    else
        %Check whether it is a logical operator
        if ~islogical(generatePlots)
            error('generatePlots function inpur must be a logical of true or false')
        end
    end

    %% Set-up

    %Set the glenoid length proportions for which the Bankart defects will
    %be created at.
    defectSizes = [0.04, 0.08, 0.12, 0.16, 0.20, 0.24];
    
    %% Calculate glenoid length
    
    % The glenoid length is used to scale the size of the defect as per
    % Klemt et al. (2019). The defect is created perpendicular to the long
    % axis of the glenoid at a proportional distance from the anterior edge
    % of the glenoid.
    
    %Create a plane to project the superior and inferior glenoid points on
    %to so that they can be measured along the same axis
    xyPlane = [0 0 0 1 0 0 0 1 0];
    
    %Project the superior and inferior glenoid points onto this plane
    projSupGlenoid = projPointOnPlane(landmarks.SupGlenoid,xyPlane);
    projInfGlenoid = projPointOnPlane(landmarks.InfGlenoid,xyPlane);
    
    %Calculate the distance between the points for glenoid length
    glenoidLength = distancePoints3d(projSupGlenoid,projInfGlenoid);
    
    %Determine anterior defect distances based on glenoid length
    defectDistances = zeros(1,length(defectSizes));
    for dd = 1:length(defectSizes)
        defectDistances(dd) = defectSizes(dd) * glenoidLength;
    end
    clear dd
    
    %% Create defect surfaces
    
    %Start the waitbar
    wbar = waitbar(0/(length(defectSizes)+1),'Creating anterior defects...');
    
    %Loop through defect sizes
    for dd = 1:length(defectSizes)
        
        %Update wait bar
        waitbar(dd/(length(defectSizes)+1),wbar,...
            ['Creating ',num2str(defectSizes(dd)*100),'% anterior bone defect...']);
        
        %Before cutting the defect, separate the top of the glenoid mesh
        %from the remainder. We can use the superior glenoid point as a
        %reference for where to cut. This just makes it easier for filling
        %holes later.
        
        %Set the cutting plane to pass through a point at the superior
        %glenoid point on the Y-axis
        cutLevelY = landmarks.SupGlenoid(2) + 3; %Set the cut level at the superior glenoid + a few mm
        snapToleranceY = mean(patchEdgeLengths(glenoidMesh.glenoidF,glenoidMesh.glenoidV))/100;
        nY = vecnormalize([0 1 0]); %Normal direction to plane
        PY = [landmarks.SupGlenoid(1) cutLevelY landmarks.SupGlenoid(3)]; %Point on plane
        
        %Cut through the glenoid at the plane (note 3rd color data output is supressed)
        [glenoidFcY,glenoidVcY,~,logicSideY,glenoidEcY] = ...
            triSurfSlice(glenoidMesh.glenoidF,glenoidMesh.glenoidV,[],PY,nY,snapToleranceY);
        
        %Visualise slice
        if generatePlots
            %Plot split planes
            cFigure; subplot(1,2,1); hold on;
            hp1 = gpatch(glenoidFcY(~logicSideY,:),glenoidVcY,'bw','none',1);
            hp2 = gpatch(glenoidFcY(logicSideY,:),glenoidVcY,'rw','none',1);
            legend([hp1 hp2],{'Surface above plane','Surface below plane'})
            axisGeom; axis manual; camlight headligth;
            colormap gjet;
            %Plot extracted surface and boundary
            subplot(1,2,2); hold on;
            gpatch(glenoidFcY(logicSideY,:),glenoidVcY,'w','none',1);
            gpatch(glenoidFcY(~logicSideY,:),glenoidVcY,'w','none',0.25);
            hp1=gpatch(glenoidEcY,glenoidVcY,'none','b',1,3);
            hp2=quiverVec(PY,nY,0.05,'k');
            legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
            axisGeom; axis manual; camlight headligth;
        end
        
        %Extract the faces we want to keep
        [glenoidKeepF,glenoidKeepV] = patchCleanUnused(glenoidFcY(logicSideY,:),glenoidVcY);
        
        %Merge vertices
        [glenoidKeepF,glenoidKeepV] = mergeVertices(glenoidKeepF,glenoidKeepV);
        
        %Fill the top part of the mesh prior to creating the defect
        glenoidTopEb = patchBoundary(glenoidKeepF,glenoidKeepV); %Get boundary edges
        indBoundary = edgeListToCurve(glenoidTopEb); %Convert boundary edges to a curve list
        indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
        angleThreshold = pi*(120/180); %threshold for self triangulation
        [glenoidKeepF,glenoidKeepV,indBoundaryTop] = ...
            triSurfSelfTriangulateBoundary(glenoidKeepF,glenoidKeepV,indBoundary,angleThreshold,1);

        %Force boundary to have a Y level aligned with the cut
        glenoidKeepV(indBoundaryTop,2) = cutLevelY;

        %Visualise the boundary on the cut surface
        if generatePlots
            cFigure; hold on;
            gpatch(glenoidKeepF,glenoidKeepV,'bw','k');
            plotV(glenoidKeepV(indBoundaryTop,:),'r-','LineWidth',2);
            camlight('headlight');
            axisGeom;
            title('Boundary of Cut Glenoid Defect')
        end

        %Create a surface that closes the glenoid defect
        %Currently uses the default 1.0 point spacing used for the glenoid
        [topF,topV] = regionTriMesh2D({glenoidKeepV(indBoundaryTop,[1 3])},1.0,0,0);
        topY = ones(length(topV),1)*mean(glenoidKeepV(indBoundaryTop,2));
        topV = [topV(:,1), topY(:,1), topV(:,2)]; %Add/set y-level converting to 3D mesh

        %Visualise new meshes
        if generatePlots
            cFigure; hold on;
            gpatch(glenoidKeepF,glenoidKeepV,'bw','k');
            gpatch(topF,topV,'gw','k');
            plotV(glenoidKeepV(indBoundaryTop,:),'r-','LineWidth',2);
            camlight('headlight');
            axisGeom;
            title('Filled Glenoid Top Surface');
        end

        %Join the two element sets
        [glenoidKeepF,glenoidKeepV,glenoidKeepC] = joinElementSets({glenoidKeepF,topF},{glenoidKeepV,topV});
        
        %Merge vertices
        [glenoidKeepF,glenoidKeepV] = mergeVertices(glenoidKeepF,glenoidKeepV);

        %Visualise the joined sets
        if generatePlots
            cFigure;
            gpatch(glenoidKeepF,glenoidKeepV,glenoidKeepC,'k');
            colormap gjet; icolorbar;
            axisGeom;
            title('Joined Surfaces for Top Filled Glenoid');
        end

        %Check boundaries to make sure there are no holes. If there is throw an
        %error as the volumetric meshing won't work.
        if ~isempty(patchBoundary(glenoidKeepF,glenoidKeepV))
            error('Holes detected in glenoid top filled mesh. Stopping here as volumetric meshing won''t work');
        end
        
        %Create the defect now by cutting along the vertical plane
        
        %Set the cutting plane to pass through a point at the specified
        %distance for the defect from the anterior glenoid point along the
        %X-axis (i.e. short axis of the glenoid), and has a normal pointing
        %along this X-axis.
        cutLevelX = landmarks.AntGlenoid(1) - defectDistances(dd); %Set the cut level
        snapToleranceX = mean(patchEdgeLengths(glenoidKeepF,glenoidKeepV))/100;
        nX = vecnormalize([1 0 0]); %Normal direction to plane
        PX = [cutLevelX landmarks.AntGlenoid(2) landmarks.AntGlenoid(3)]; %Point on plane
        
        %Cut through the glenoid at the plane (note 3rd color data output is supressed)
        [glenoidFcX,glenoidVcX,~,logicSideX,glenoidEcX] = ...
            triSurfSlice(glenoidKeepF,glenoidKeepV,[],PX,nX,snapToleranceX);
        
        %Visualise slice
        if generatePlots
            %Plot split planes
            cFigure; subplot(1,2,1); hold on;
            hp1 = gpatch(glenoidFcX(~logicSideX,:),glenoidVcX,'bw','none',1);
            hp2 = gpatch(glenoidFcX(logicSideX,:),glenoidVcX,'rw','none',1);
            legend([hp1 hp2],{'Surface above plane','Surface below plane'})
            axisGeom; axis manual; camlight headligth;
            colormap gjet;
            %Plot extracted surface and boundary
            subplot(1,2,2); hold on;
            gpatch(glenoidFcX(logicSideX,:),glenoidVcX,'w','none',1);
            gpatch(glenoidFcX(~logicSideX,:),glenoidVcX,'w','none',0.25);
            hp1=gpatch(glenoidEcX,glenoidVcX,'none','b',1,3);
            hp2=quiverVec(PX,nX,0.05,'k');
            legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
            axisGeom; axis manual; camlight headligth;
        end
        
        %Extract the faces we want to keep
        [glenoidDefectF,glenoidDefectV] = patchCleanUnused(glenoidFcX(logicSideX,:),glenoidVcX);

        %Merge vertices
        [glenoidF,glenoidV] = mergeVertices(glenoidDefectF,glenoidDefectV);

        %Visualise glenoid defect as a whole
        if generatePlots
           cFigure
           gpatch(glenoidF,glenoidV,'kw','none');
           axisGeom; camlight headlight
           title([num2str(defectSizes(dd)*100),'% Anterior Glenoid Defect'])
        end
        
        %Self triangulate the potentially jagged edge of the cut
        glenoidDefectEb = patchBoundary(glenoidF,glenoidV); %Get boundary edges
        
        %Depending on the cut level and the shape of the glenoid, more than
        %one boundary can be create. Test this out by grouping boundaries,
        %and looking for whether there is more than one.
        optionStruct.outputType='label';
        G = tesgroup(glenoidDefectEb,optionStruct);
        
        %Check if there is more than one boundary
        if max(G) > 1
            
            %Indicates that there is more than one boundary
            %Idenfity how many boundaries there are
            numEb = max(G);
            
            %Loop through boundaries and fill
            for bb = 1:numEb
                
                %Get the current boundary
                glenoidDefectEb_curr = glenoidDefectEb(G == bb,:);
                indBoundary = edgeListToCurve(glenoidDefectEb_curr); %Convert boundary edges to a curve list
                indBoundaryDefect = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
% % %                 angleThreshold = pi*(120/180); %threshold for self triangulation
% % %                 [glenoidF,glenoidV,indBoundaryDefect] = ...
% % %                     triSurfSelfTriangulateBoundary(glenoidF,glenoidV,indBoundary,angleThreshold,1);

                %Force boundary to have an X level aligned with the cut
                glenoidV(indBoundaryDefect,1) = cutLevelX;
                
% % %                 %The hard edge of the cut at 90 degrees causes convergence
% % %                 %issues when simulating the humeral head surface over the
% % %                 %glenoid. A potential solution here is to create a chafered
% % %                 %edge. To do this we need to both contract the size of the
% % %                 %boundary and shift it along the x-axis.
% % %                 
% % %                 %Push to boundary 1mm away from the edge of the cut
% % %                 defectPoints = glenoidV(indBoundaryDefect,:);
% % %                 defectPoints(:,1) = cutLevelX + 1.00;
% % %                 
% % %                 %Contract the size of the defect by 1mm
% % %                 %Set up polygon
% % %                 pg = [defectPoints(:,2),defectPoints(:,3)];
% % %                 %Contract polygon
% % %                 %Uses positive number for some reason???
% % %                 %These values need to be more rigorously tested as they are
% % %                 %somewhat arbitrary to work for the current examp;e
% % %                 pgIn = expandPolygon(pg, 0.25, 'cleanuploops', true);
% % %                 %Set back into the defect points structure. Throw error if
% % %                 %mismatching length
% % %                 if length(pgIn{1,1}) ~= length(defectPoints)
% % %                     error('Mismatching length between defect points and contracted polygon.')
% % %                 end
% % %                 defectPointsContracted = [defectPoints(:,1),pgIn{1,1}(:,1),pgIn{1,1}(:,2)];
                
                %Visualise the boundary on the cut surface
                if generatePlots
                    cFigure; hold on;
                    gpatch(glenoidF,glenoidV,'bw','k');
                    plotV(glenoidV(indBoundaryDefect,:),'r-','LineWidth',2);
% % %                     plotV(defectPointsContracted,'r-','LineWidth',2);
                    camlight('headlight');
                    axisGeom;
                    title('Boundary of Cut Glenoid Defect')
                end

                %Create a surface that closes the glenoid defect
                %Currently uses the default 1.0 point spacing used for the glenoid
                [defectF,defectV] = regionTriMesh2D({glenoidV(indBoundaryDefect,[2 3])},1.0,0,0);
% % %                 [defectF,defectV] = regionTriMesh2D({defectPointsContracted(:,[2 3])},1.0,0,0);
                defectX = ones(length(defectV),1)*mean(glenoidV(indBoundaryDefect,1));
% % %                 defectX = ones(length(defectV),1)*mean(defectPointsContracted(:,1));
                defectV = [defectX(:,1), defectV(:,1), defectV(:,2)]; %Add/set x-level converting to 3D mesh

                %Visualise new meshes
                if generatePlots
                    cFigure; hold on;
                    gpatch(glenoidF,glenoidV,'bw','k');
                    gpatch(defectF,defectV,'gw','k');
                    plotV(glenoidV(indBoundaryDefect,:),'r-','LineWidth',2);
                    camlight('headlight');
                    axisGeom;
                    title('Filled Glenoid Defect Surface');
                end
                
                %Join the two element sets
                [glenoidF,glenoidV,glenoidC] = joinElementSets({glenoidF,defectF},{glenoidV,defectV});
                
                %Visualise the boundary on the cut surface
                if generatePlots
                    cFigure; hold on;
                    gpatch(glenoidF,glenoidV,'bw','k');
                    plotV(glenoidV(indBoundary,:),'r-','LineWidth',2);
                    camlight('headlight');
                    axisGeom;
                    title('Boundary of Cut Glenoid Defect')
                end
                
% % %                 %Sweep loft between the original defect and contracted
% % %                 %curve
% % %                 numStepsCurve = 3; %Number of steps for the curve
% % %                 p1 = mean(glenoidV(indBoundary,:),1); %First point
% % %                 n1 = [1,0,0]; %First direction vector
% % %                 p2 = mean(defectPointsContracted,1); %End point
% % %                 n2 = [1,0,0]; %End direction vector
% % %                 csapsSmoothPar = 1; %Cubic smoothening spline smoothening parameter
% % %                 f = 0.05; %Extent of tangential nature to boundary curves, surface will remain approximately orthogonal to input curves for f*distance between curves
% % %                 [Vg] = sweepCurveSmooth(p1,p2,n1,n2,numStepsCurve,csapsSmoothPar,f);
% % %                 [F,V,C] = sweepLoft(glenoidV(indBoundaryDefect,:),defectPointsContracted,n1,n2,Vg);
% % %                 
% % %                 %Convert quad loft to tri
% % %                 [Ft,Vt,Ct] = quad2tri(F,V,'f',C);
% % %                 
% % %                 %Visualise
% % %                 if generatePlots
% % %                     cFigure; hold on;
% % %                     gpatch(glenoidF,glenoidV,'bw','k');
% % %                     plotV(glenoidV(indBoundary,:),'r-','LineWidth',2);
% % %                     gpatch(Ft,Vt,'gw','k');
% % %                     camlight('headlight');
% % %                     axisGeom;
% % %                     title('Boundary of Cut Glenoid Defect')
% % %                 end
% % %                 
% % %                 %Join the element sets
% % %                 [glenoidF,glenoidV,glenoidC] = joinElementSets({glenoidF,Ft},{glenoidV,Vt});
                                
            end
            clear bb

        else
            
            %Only one boundary is present        
            indBoundary = edgeListToCurve(glenoidDefectEb); %Convert boundary edges to a curve list
            indBoundaryDefect = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
% % %             angleThreshold = pi*(120/180); %threshold for self triangulation
% % %             [glenoidF,glenoidV,indBoundaryDefect] = ...
% % %                 triSurfSelfTriangulateBoundary(glenoidF,glenoidV,indBoundary,angleThreshold,1);
            
            %Force boundary to have an X level aligned with the cut
            glenoidV(indBoundaryDefect,1) = cutLevelX;

            %Visualise the boundary on the cut surface
            if generatePlots
                cFigure; hold on;
                gpatch(glenoidF,glenoidV,'bw','k');
                plotV(glenoidV(indBoundaryDefect,:),'r-','LineWidth',2);
                camlight('headlight');
                axisGeom;
                title('Boundary of Cut Glenoid Defect')
            end

            %Create a surface that closes the glenoid defect
            %Currently uses the default 1.0 point spacing used for the glenoid
            [defectF,defectV] = regionTriMesh2D({glenoidV(indBoundaryDefect,[2 3])},1.0,0,0);
            defectX = ones(length(defectV),1)*mean(glenoidV(indBoundaryDefect,1));
            defectV = [defectX(:,1), defectV(:,1), defectV(:,2)]; %Add/set x-level converting to 3D mesh

            %Visualise new meshes
            if generatePlots
                cFigure; hold on;
                gpatch(glenoidF,glenoidV,'bw','k');
                gpatch(defectF,defectV,'gw','k');
                plotV(glenoidV(indBoundaryDefect,:),'r-','LineWidth',2);
                camlight('headlight');
                axisGeom;
                title('Filled Glenoid Defect Surface');
            end

            %Join the two element sets
            [glenoidF,glenoidV,glenoidC] = joinElementSets({glenoidF,defectF},{glenoidV,defectV});
            
        end

        %Merge vertices for fileld glenoid
        [glenoidF,glenoidV] = mergeVertices(glenoidF,glenoidV);

        %Visualise the joined sets
        if generatePlots
            cFigure;
            gpatch(glenoidF,glenoidV,glenoidC,'k');
            colormap gjet; icolorbar;
            axisGeom;
            title('Joined Surfaces for Cut and Filled Glenoid');
        end

        %Check boundaries to make sure there are no holes. If there is throw an
        %error as the volumetric meshing won't work.
        if ~isempty(patchBoundary(glenoidF,glenoidV))
            error('Holes detected in glenoid defect mesh. Stopping here as volumetric meshing won''t work');
        end
        
        %Remesh the glenoids to the desired point spacing
        if ~isempty(pointSpacing)
            [glenoidFb,glenoidVb] = triRemeshLabel(glenoidF,glenoidV,pointSpacing);
        end

        %Output the created meshes
        varName = ['anteriorDefect_',num2str(defectSizes(dd)*100),'per'];
        glenoidDefectMesh.(char(varName)).glenoidF = glenoidFb;
        glenoidDefectMesh.(char(varName)).glenoidV = glenoidVb;
        
    end
    clear dd
    
    %% Finish process
    
    %Update waitbar
    waitbar((length(defectSizes)+1)/(length(defectSizes)+1),wbar,...
            'Finished creating anterior bone defects!');
    
    %Close waitbar
    close(wbar)

%% %%%%% ----- End of createAnteriorBankartDefects.m ----- %%%%%