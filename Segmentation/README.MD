# Segmentation & Surface Meshing Processes

The following sections outline the details and processes used for segmenting and eatracting surface meshes from the Visible Human Dataset. Thoughout the folders included in this directory (i.e. Mimics, 3matic and Exported) will be referred to.

## Segmentation in Mimics

Segmentation of the CT data was performed in Materialise Mimics Research v21.0. The major tool used was the threshold feature to create masks of the bony elements of the scans. Mimics splitting tool was used to extract indiviual parts. The various filling tools were used to fill out the masks and create subsequent complete surface parts of the relevant bodies. Manual refinement of the masks was performed using the other editing tools (e.g. threshold painting) to finalise the parts.

***TODO:***

- Add details of bodies extracted and process to do so.

## Surface Meshing in 3matic

Meshing of the parts was a slightly more complex process than segmentation, and was performed across a number of steps. These will be outlined here. Materialise 3matic Research v13.0 was used to export surface meshes that were post processed with the [GIBBON toolbox](https://www.gibboncode.org/) in Matlab and will be mentioned throughout. 

### Meshing Steps

The following section outlines the steps for creating the meshes used for analysis:

- The masks were originally converted to 'Parts' in Mimics and copied across to the 3matic project workspace.
- The 3matic 'Smooth' feature was used with a smooth factor of 0.7 to smooth the surfaces
- The 3matic 'Uniform Remesh' feature was used to get uniform triangle sizes of 1.5 on both the humerus and scapula surfaces. Sharp angles (60.0) were preserved, but not really expected. The triangle edge size was based off that used by Moroder et al. (2019) who used a glenoid part in a similar manner
- Following these steps a new 3matic project file was created to rotate the parts to a *'glenoid coordinate system'*.
- A 3D plane was fit to the glenoid by marking the triangles on the glenoid face surface and fitting a plane to this. This plane will serve as the XY plane of the glenoid coordinate system. 
- A point was added to the deepest point of the glenoid which will serve as the glenoid coordinate systems origin. The centre of the plane fitted to the glenoid surface assisted in identifying this point.
- The supra and infraglenoid tubercles were marked with points, which will serve as the long axis of the glenoid (i.e. Y-axis).
- A plane was fit that passed through the two tubercle points and was perpendicular to the glenoid surface plane. This plane will serve as the YZ plane in the glenoid coordinate system.
- Original world coordinate system planes were created for the XY, XZ and YZ planes in the original coordinate system. This was done by creating planes that ran through three points on thes planes (e.g. using relevant 0 and 1 coordinates).
- The 'Plane to Plane Align' function of 3matic was used to make the glenoid surface plane coincident with the world XY plane, thus making the vector pointing away from the glenoid surface the Z-axis.
- The same 'Plane to Plane Align' feature was then repeated to align the plane fit to the tubercles to the world YZ plane, making the line between the tubercles the Y-axis of the glenoid coordinate system.
- The 'Translate' function of 3matic was then used to shift all of the objects and make the deep glenoid point the origin (i.e. 0,0,0).
- Following all of this, the X-axis essentially also corresponds with the short axis of the glenoid. 
- After running createGlenoidCartilage.m in MATLAB, a .stl file called 'GlenoidCarilage_GlenoidCoordinateSystem.stl' is created. This is imported back into the GlenoidCoordinateSystem 3matic file as a new structure (Import Part; set Fix Normals to be checked and scale coefficient in mm).
- Next we fit a sphere to the humeral head to approximate the humeral head size and position. This was achieved by marking triangles around the humeral head and fitting a sphere analytical primitive to the marked triangles.
- We used the marked triangles in a similar manner to above in creating articular cartilage on the humeral head, by extracting these triangles as a surface and smoothing the edges (filter distance = 0.3 for harp triangles). This STL surface was then exported and used in the createHumeralHeadCartilage.m function (see Code > Main folder). 
- Like earlier this humeral head cartilage can be imported into the curretn 3matic file (same parameters as before).
- The humerus is clearly rotated and offset from the scapula in it's segmented form (i.e. not a neutral glenohumeral position). We create a new 3matic project file that will have an aligned neutral humerus position.
- We created the coordinate system axes for the scapula and humerus bodies as per Wu et al. (2005). In the absence of the elbow epicondyles on the humerus, we had to approximate certain aspects of the humeral coordinate system. Specifically, the longitudinal axes was a line fit to the centre of the humeral shaft; and a lateal point on the greater tubercle to indicate the lateral directional axis. After defining the coordinate systems, the appropriate rotations were applied to the humerus body to align all three of the axes with the scapula coordiante system. To do this, three planes were crated for each coordinate system and these were subsequently aligned with one another.
- After alignment was completed we saved a new 3matic project file to create the base system for the FE analysis.
- We fit a sphere to the inside aspect of the created humeral head cartilage to create an *'idealised huemral head'* that would serve as a rigid body for the simulations. This was created by convering the analytical primitive sphere to a part.
- The next step is to create the glenoid surface cartilage. We carefully marked the triangles that were on the glenoid surface as a means to create an initial rough outline edge of the glenoid surface. We create a copy of the scapula so that we can smooth the border of these marked triangles without affecting the original mesh. Smoothing the edges of the marked triangles is achieved through the 'Smooth Marking Border' in 3matic. After this we now have a smooth edge surface of the glenoid face that we can extract from the scapula object to its own object (under the Surface List tree in 3matic). This surface can now be exported as an STL. Prior to processing, it is beneficial to iteratively mark the face triangles and smooth the border multiple times to ensure a smooth part edge. Lastly, filtering sharp triangles (filter distance = 0.4) produced a smooth edge by increasing triangle size around the edges evenly. Further processing in the next stop is done with the GIBBON toolbox using the createGlenoidCartilage.m function (see Code folder) (***TODO: check if this holds as the process for glenoid cartilage***).
- The surfaces of the humerus and scapula are also exported at this point as they are needed to create the glenoid cartilage mesh (if using the constant thickness method).
- We recreated a smoother version of the humeral head cartilage on this sphere using a slightly different process to above. A circle was fit that approximated the border of the humeral cartilage and this was used to cut the idealised humeral head to extract the relevant surface for articular cartilage. The cutting entity was defnied as a plane that was parallel to the glenoid cartilage arc created. The STL of this is exported in a similar fashion as before and used in the createIdealisedHumeralHeadCartilage.m function (see Code folder). The STL made for the idealised cartilage was then imported back in to 3matic (***TODO: check if this holds as the process for humeral head cartilage - matter of fact check the whole process here to make sure cartilage process is correct...***).
- We created a sphere at the same centre point as the iealised humeral head, but with a larger radius relevant to the curvature radius of the cartilage. The radius of this 'cartilage sphere' was increased by the same proportions from humeral head:humeral cartilage radius of curvature presented in Walia et al. (2013; 2015). We translated this new sphere back along the normal vector of the plane we used to cut the idealised humeral head sphere into its cartilage section by a distance that was the difference of the two sphere radii, so that the sphere edges coincided at a point directly away form the normal of the cutting plane. We then shifted the sphere back along this normal direction by 2.03mm to replicate the minimum thickness of the humeral cartilage also specified in Walia et al. (2013; 2015). Lastly, we cut this larger sphere using the same plane we used to cut the humeral head sphere into its cartilage surface. From this, we not have two cut spheres, one slightly smaller representing the base surface of the cartilage, and the other larger one representing the outer surface of the cartilage. These could then be used in the 'createCartilageSurfaces' Matlab code (see Code > Main) to create the humeral head cartilage.
- The created cartilages were imported back into the 'BaseSystem' 3matic file. As part of this we inverted or fixed the face normals to be pointing out. 
- We then offset the humeral head aspects away from the glenoid by 5mm. This potentially presented an 'unnatural' position of the joint, but provided a consistent starting position with whch simulations could be run.

### 3matic File System

There are a number of different 3matic project files included in the firectory. These are explained here:

**1_BaselineMeshing:** This contains the baseline smoothed and remeshed humerus and scapula in their original coordinate system corresponding to the CT scans.

**2_GlenoidCoordinateSystem:** This contains the humerus and scapula meshes rotated to the *'glenoid coordinate system'* whereby the deep glenoid point is the origin, the z-axis points directly away from the glenoid surface plane, and the x- and y-axis represent the horizontal (i.e. 3 o'clock to 9 o'clock) and vertical axes (i.e. 12 o'clock to 6 o'clock; or supra to infraglenoid tubercle) of the glenoid.

**3_AlignedHumerus:** This builds on the glenoid coordinate system file and rotates/translates the humerus into what we would consider a neutral glenohumeral position (i.e. 0 degrees of axial rotation and elevation).

**4_BaseSystem:** This builds on the glenoid coordinate system and aligned humerus files by creating the base system for the FE analysis which includes the idealised humeral head sphere and it's cartilage surface, along with the simplified detached glenoid section.


## TODO:
- Complete the list of different project files and what they contain (i.e. baseline file, rotated file etc.)
- Add any remaining references

## References

Moroder et al. (2019). Challenging the current concept of critical glenoid bone loss in shoulder instability. Does the size measurement really tell all? *Am J Sports Med*, 47: 688-694.

Walia et al. (2013). Theoretical model of the effect of combined glenohumeral bone defects on anterior shoulder instability: A finite element approach. *J Orthop Res*, 31: 601-607.

Walia et al. (2015). Influence of combined Hill-Sachs and bony Bankart defects on range of motion in anterior instability in a finite element model. *Arthroscopy*, 31: 2119-2127.

Wu et al. (2005). ISB recommendation on definitions of joint coordinate systems of various joints for the reporting of human joint motion - Part II: Shoulder, elbow, wrist and hand. *J Biomech*, 38: 981-992.