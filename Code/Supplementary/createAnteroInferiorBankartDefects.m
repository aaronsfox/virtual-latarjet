function [glenoidDefectMesh] = createAnteroInferiorBankartDefects(glenoidMesh,landmarks,generatePlots)
    
    %% TODO: consider inlcluding Latarjet process with this? Or separate

    %% This function serves to import in and create the base surface system for
    %  running the FEA analysis of the humeral head against the glenoid.
    %
    %  Inputs:
    %
    %   glenoidMesh         structure of glenoid mesh containing faces and
    %                       vertices as 'glenoidF' and 'glenoidV'
    %   landmarks           structure containing glenoid landmarks 
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
    % References
    %
    % Itoi et al. (2000). Effect of a glenoid defect on anteroinferior
    % stability of the shoulder after Bankart repair. J Bone Joint Surg,
    % 82-A: 35-46.
    %
    % Walia et al. (2013). Theoretical model of the effect of combined
    % glenohumeral bone defects on anterior shoulder instability: A finite
    % element approach. J Orthop Res, 31: 601-607.
    %
    % Walia et al. (2015). Influence of combined Hill-Sachs and bony
    % Bankart defects on range of motion in anterior instability of the
    % shoulder in a finite element model. Arthroscopy, 31: 2119-2127.
    
    %% Check inputs

    %Participant directory
    if nargin < 2
        error('Structures containing the glenoid mesh as well as landmarks are required')
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

    %% Set-up

    %Set the defect sizes to create. These values are relative to the
    %radius of a circle fitted to the glenoid length
    defectSizes = [0.25, 0.50, 0.75, 1.00];
    
    %% Calculate glenoid length
    
    % The glenoid length is used to scale the size of the defect as per
    % Itoi et al (2000) and Walia et al. (2013, 2015). A circle is fitted
    % to match glenoid length and the radius of this is used to scale the
    % anteroinferior defects.
    
    %Create a plane to project the superior and inferior glenoid points on
    %to so that they can be measured along the same axis
    xyPlane = [0 0 0 1 0 0 0 1 0];
    
    %Project the superior and inferior glenoid points onto this plane
    projSupGlenoid = projPointOnPlane(landmarks.SupGlenoid,xyPlane);
    projInfGlenoid = projPointOnPlane(landmarks.InfGlenoid,xyPlane);
    
    %Calculate the distance between the points for glenoid length
    glenoidLength = distancePoints3d(projSupGlenoid,projInfGlenoid);
    
    %Determine the size of the defects to be created based on the radius
    %(i.e. half) of glenoid length
    defectDistances = zeros(1,length(defectSizes));
    for dd = 1:length(defectSizes)
        defectDistances(dd) = defectSizes(dd) * (glenoidLength/2);
    end
    clear dd
    
    %% Create defect surfaces
    
    %Start the waitbar
    wbar = waitbar(0/(length(defectSizes)+1),'Creating anteroinferior defects...');
    
    %Loop through defect sizes
    for dd = 1:length(defectSizes)
        
        %Update wait bar
        waitbar(dd/(length(defectSizes)+1),wbar,...
            ['Creating ',num2str(defectSizes(dd)*100),'% glenoid radius anteroinferior bone defect...']);
        
        %To create the defect cutting plane we must generate a plane that
        %is at a 45 degree angle to the glenoid at the glenoid length
        %radius with the defect distance subtracted. Given the glenoid is
        %oriented en face to the XY plane, we can simply use these axes as
        %the normal direction for a 45 degree angled plane --- and take the
        %value from the deep glenoid origin in this direction to define the
        %plane.
        
        %Identify the mid point of the superior and inferior glenoid points
        %along the y axis. This mid point might be slightly different to
        %the deep glenoid point, so is necessary for fitting the circle
        %described in Itoi et al. and Walia et al. This also uses some
        %other relevant points to define this in 3d.
        vertMidPt = [0 (projSupGlenoid(2)+projInfGlenoid(2))/2 landmarks.SupGlenoid(3)];
        
        %Define point to cut through based on the glenoid length and
        %current defect size
        
        %Determine horizontal (X) and vertical (-Y) distance to the point
        %that needs to be cut (i.e. along the glenoid radius - defect
        %distance)
        horzDist = cosd(45) * (glenoidLength/2 - defectDistances(dd));
        vertDist = sind(45) * (glenoidLength/2 - defectDistances(dd)) * -1;
        
        %Determine point for angled plane to pass through
        ptTransform = createTranslation3d(horzDist, vertDist, 0);
        P = transformPoint3d(vertMidPt, ptTransform);
        
        %Set the normal for the cutting plane (i.e. 45 degrees to XY plane)
        n = vecnormalize([1 -1 0]);

        %Set snap tolerance for the cutting function
        snapTolerance = mean(patchEdgeLengths(glenoidMesh.glenoidF,glenoidMesh.glenoidV))/100;
        
        %Cut through the glenoid at the plane (note 3rd color data output is supressed)
        [glenoidFc,glenoidVc,~,logicSide,glenoidEc] = ...
            triSurfSlice(glenoidMesh.glenoidF,glenoidMesh.glenoidV,[],P,n,snapTolerance);
        
        %Visualise slice
        if generatePlots
            %Plot split planes
            cFigure; subplot(1,2,1); hold on;
            hp1 = gpatch(glenoidFc(~logicSide,:),glenoidVc,'bw','none',1);
            hp2 = gpatch(glenoidFc(logicSide,:),glenoidVc,'rw','none',1);
            legend([hp1 hp2],{'Surface above plane','Surface below plane'})
            axisGeom; axis manual; camlight headligth;
            colormap gjet;
            %Plot extracted surface and boundary
            subplot(1,2,2); hold on;
            gpatch(glenoidFc(logicSide,:),glenoidVc,'w','none',1);
            gpatch(glenoidFc(~logicSide,:),glenoidVc,'w','none',0.25);
            hp1=gpatch(glenoidEc,glenoidVc,'none','b',1,3);
            hp2=quiverVec(P,n,0.05,'k');
            legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
            axisGeom; axis manual; camlight headligth;
        end
        
        %Extract the faces we want to keep
        [glenoidKeepF,glenoidKeepV] = patchCleanUnused(glenoidFc(logicSide,:),glenoidVc);
        
        %Merge vertices
        [glenoidF,glenoidV] = mergeVertices(glenoidKeepF,glenoidKeepV);

        %Visualise glenoid defect as a whole
        if generatePlots
           cFigure
           gpatch(glenoidF,glenoidV,'kw','none');
           axisGeom; camlight headlight
           title([num2str(defectSizes(dd)*100),'% glenoid radius anteroinferior bone defect'])
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
            %This shouldn't really happen with the anteroinferior defects,
            %maybe at the largest if the participant has a prominent
            %coracoid process...
            %Idenfity how many boundaries there are
            numEb = max(G);
            
            %Rotate the glenoid vertices so that the cut plane is aligned
            %with the xy plane
            
            %Create the transform plane
            cutPlane = createPlane(P,n);
            glenoidTransform = createBasisTransform3d(xyPlane,cutPlane);
            
            %Transform glenoid surface
            for pp = 1:length(glenoidV)
                glenoidV(pp,:) = transformPoint3d(glenoidV(pp,:),glenoidTransform);    
            end
            clear pp
            
            %Loop through boundaries and fill
            for bb = 1:numEb
                
                %Get the current edge boundary
                glenoidDefectEb_curr = glenoidDefectEb(G == bb,:);
                indBoundary = edgeListToCurve(glenoidDefectEb_curr); %Convert boundary edges to curve list
                indBoundaryDefect = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop

                %Force boundary to have a Z-level of 0 now that it's been
                %rotated to the XY plane
                glenoidV(indBoundaryDefect,3) = 0;

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
                [defectF,defectV] = regionTriMesh2D({glenoidV(indBoundaryDefect,[1 2])},1.0,0,0);
                defectZ = ones(length(defectV),1)*mean(glenoidV(indBoundaryDefect,3));
                defectV = [defectV(:,1), defectV(:,2), defectZ(:,1)]; %Add/set z-level converting to 3D mesh

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
                
            end
            clear bb
            
            %Rotate the joined element sets back to the original position
            reverseTransform = createBasisTransform3d(cutPlane,xyPlane);
            for pp = 1:length(glenoidV)
                glenoidV(pp,:) = transformPoint3d(glenoidV(pp,:),reverseTransform);    
            end
            clear pp

        else
            
            %Only one boundary is present  
            %This should typically happen with the anteroinferior defects
            
            %Rotate the glenoid vertices so that the cut plane is aligned
            %with the xy plane
            
            %Create the transform plane
            cutPlane = createPlane(P,n);
            glenoidTransform = createBasisTransform3d(xyPlane,cutPlane);
            
            %Transform glenoid surface
            for pp = 1:length(glenoidV)
                glenoidV(pp,:) = transformPoint3d(glenoidV(pp,:),glenoidTransform);    
            end
            clear pp
            
            %Get edge boundery
            glenoidDefectEb = patchBoundary(glenoidF,glenoidV); %Get boundary edges
            indBoundary = edgeListToCurve(glenoidDefectEb); %Convert boundary edges to a curve list
            indBoundaryDefect = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop

            %Force boundary to have a Z-level of 0 now that it's been
            %rotated to the XY plane
            glenoidV(indBoundaryDefect,3) = 0;

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
            [defectF,defectV] = regionTriMesh2D({glenoidV(indBoundaryDefect,[1 2])},1.0,0,0);
            defectZ = ones(length(defectV),1)*mean(glenoidV(indBoundaryDefect,3));
            defectV = [defectV(:,1), defectV(:,2), defectZ(:,1)]; %Add/set z-level converting to 3D mesh

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
            
            %Rotate the joined element sets back to the original position
            reverseTransform = createBasisTransform3d(cutPlane,xyPlane);
            for pp = 1:length(glenoidV)
                glenoidV(pp,:) = transformPoint3d(glenoidV(pp,:),reverseTransform);    
            end
            clear pp
            
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
        
        %Output the created meshes
        varName = ['anteroInferiorDefect_',num2str(defectSizes(dd)*100),'per'];
        glenoidDefectMesh.(char(varName)).glenoidF = glenoidF;
        glenoidDefectMesh.(char(varName)).glenoidV = glenoidV;
        
    end
    clear dd
    
    %% Finish process
    
    %Update waitbar
    waitbar((length(defectSizes)+1)/(length(defectSizes)+1),wbar,...
            'Finished creating anteroinferior bone defects!');
    
    %Close waitbar
    close(wbar)

%% %%%%% ----- End of createAnteroInferiorBankartDefects.m ----- %%%%%