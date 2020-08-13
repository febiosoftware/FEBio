# Roadmap

This roadmap briefly discusses some planned milestones for the FEBio project. 

Please check out our [contribution guidelines](CONTRIBUTING.md) and [code of conduct](CODE_OF_CONDUCT.md) for information on how to contribute to the FEBio project.

## Main Milestones

The main focus at this point is the implementation of the main aims of the NIH RO1 grant that funds FEBio development. The estimated completion date for these milestones is summer of 2023, although some functionality may be available sooner.

### Formulate and implement a computationally efficient damage and fatigue failure framework for fibrous tissues. 

This aim focuses on the development and implementation of a novel theoretical framework for modeling damage and fatigue failure mechanics. This work will extend the remodeling and growth framework already implemented in FEBio to include damage and failure of soft tissues. 

### Expand the multiphysics capabilities of FEBio to include solute transport, reactions involving solutes, and tissue growth and remodeling in fluid domains at their interfaces. 

In this aim, the multiphisics capabilities of FEBio will be extended to include various transport and reaction modeling capabilites at solid-fluid boundaries. 

### Integrate the use of image data through the entire simulation pipeline, from model setup to model validation. 

This aim proposes to include image data throughout the entire modeling pipeline, from pre-processing to post-processing. Most of this development will actually occur in FEBio Studio, but the focus in FEBio is the representation of heterogeneous model parameters (e.g. material parameters) that are extracted from image data. 

## Long term milestones

The following topics are currently under consideration as potential milestones, depending on user interest and feasibility. 

### topology optimization

FEBio already has optimization algorithms for optimizing model parameters while preserving the geometry. In this milestone we would look into extending this framework to allow for topology changes. 

### Integration with the Julia programming language. 

The Julia programming language appears to be an interesting potential path for extending FEBio's capabilities through scripting, or as a mechanism for integrating FEBio with other tools. 

