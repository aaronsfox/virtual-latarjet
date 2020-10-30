# Segmentation & Surface Meshing Processes

The following sections outline the details and processes used for segmenting and extracting surface meshes from the [NMDID Dataset](https://nmdid.unm.edu/). Each participant has their own folder which includes the Mimics segmentation file, the 3matic project file, and the extracted surface meshes of the humeru and scapula.

## Segmentation in Mimics

Segmentation of the CT data was performed in Materialise Mimics Research v21.0. The major tool used was the threshold feature to create masks of the bony elements of the scans. Mimics splitting tool was used to extract indiviual parts. The various filling tools were used to fill out the masks and create subsequent complete surface parts of the relevant bodies. Manual refinement of the masks was performed using the other editing tools (e.g. threshold painting) to finalise the parts. The file was saved in the 'Segmentation' folder in a folder with the participant ID as 'id*_Segmentation.mcs'.

The masks of the scapula and humerus were converted to parts and copied across to a 3matic project file, where some further processing (see Surface Meshing notes below) were undertaken before exporting as an ASCII STL file.

## Surface Meshing in 3matic

Some processing of the meshes prior to exporting was performed in 3matic. This included some surface smoothing as well as marking of landmarks for processing later on. Materialise 3matic Research v13.0 was used to edit and export surface meshes. 

### Meshing Steps

The following section outlines the steps for creating the meshes used for analysis:

- The 3matic 'Smooth' feature was used with a smooth factor of 0.7 to smooth the scapula and humerus surfaces.
- The 3matic 'Adaptive Remesh' feature was used to remesh the triangles of the surfaces. A target edge length of 0.5mm and maximum edge length of 1.0mm were used while maintaining surface contours. This resulted in quite a fine mesh but this will be re-assessed later with mesh refinement processes.
- A sphere was fit to the humeral head ('HeadSphere') by marking selected triangles with the brush mark and fitting an analytical sphere to these. The triangles were copiedand retained in a new surface part to be exported later.
- A series of landmarks were marked according to the definitions provided by Wu et al. (2005).
- A 'HHC' point was created at the centre of the aforementioned sphere.
- The epicondyle points of the humerus were marked and labelled as 'EL' and 'EM'.
- The angulus acromialis ('AA'), trigonum scapulae ('TS'), angulus inferior ('AI'), the most dorsal point of the acromioclavicular joint ('AC'), and most ventral point on the coracoid process ('PC') were marked on the scapula.
- The supraglenoid ('SGT') and infraglenoid ('IGT') tubercles of the scapula were marked.
- The superior ('SupGlenoid'), inferior ('InfGlenoid'), anterior ('AntGlenoid') and posterior ('PostGlenoid') points of the glenoid rim were estimated and marked (i.e. widest points).
- A plane was fit to the glenoid surface ('GlenoidPlane') by marking triangles on the glenoid face. These triangles were copied and retained as a separate part.
- This plane was used to identify and create a point at the central deepest point of the glenoid surface ('DeepGlenoidPoint').
- Prior to exporting, the humerus and scapula surfaces were remeshed using the 3matic 'Uniform Remesh' function (target triangle edge length = 0.5, preserve sharp edges = true, sharp edge angle = 60). This is quite a fine mesh, however we will adjust and potentially adapt this with a later mesh refinement procedure.
- As much of the coracoid section prior to hitting scapula/glenboid was marked, and a cylinder fit to this ('CoracoidCylinder'). The axis of this cylinder was used as the normal of a line from the coracoid tip (i.e. end of cylinder). We can create a line starting at the coracoid tip and prjecting along the cylinder axis by a specified distance. The distance was based off two factors . Firstly, guided by the 22 to 27mm recommendation of Bhatia et al. (2014). Secondly, only going to a point that still is onyl coracoid, rather than projecting into the scapula. This line was created and saved as a landmark ('GraftLine'). A plane was then fit using the origin as the end of this line, and the direction of this line as the normal ('GraftPlane').
- To create the plane to shave the inferior surface of the coracoid for the standard graft, an estimate of the inferior surface of the coracoid was marked and a plane was fit to this ('InferiorPlaneFit'). This plane may be OK, but to maximise the size of the graft we use a line along the normal of this plane to a point on the coracoid that creates a flat surface, while maximising the size of the graft ('InferiorShaveLine'). A plane was then fit at the end of this line, parallel to the original fit plane ('InferiorShavePlane').

***TODO: figure out how to create planes in the right spot and in the right orientation...***

- The same process as above was followed but just using te anterior portion of the coracoid to create the landmarks/planes for te congruent arc technique (i.e. 'AnteriorPlaneFit', 'AnteriorShaveLine', 'AnteriorShavePlane').
- The humerus ('Humerus_r') and scapula ('Scapula_r'), along with the extracted triangles of the humeral head surface ('HumerusArticulating_r') and glenoid surface ('GlenoidSurface_r'), were exported as ASCII .stl files.
- All landmarks and objects were exported as 'Curves and Primitives' to separate .xml files.
- All files were saved in the participant's ID folder under the 'Segmentation' folder.
- These files were used as inputs in a subsequent Matlab code (see: ***TODO: add .m filename***).

## References

Bhatia et al. (2014). The outcomes and surgical techniques of the Latarjet procedure. *Arthroscopy*, 30: 227-235. 

Wu et al. (2005). ISB recommendation on definitions of joint coordinate systems of various joints for the reporting of human joint motion - Part II: Shoulder, elbow, wrist and hand. *J Biomech*, 38: 981-992.