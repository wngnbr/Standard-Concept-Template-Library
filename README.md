# Standard Concept Template Library

This repository contains official implementations of the Standard Concept Template Library (STL-C) associated with the NeurIPS 2024 Datasets and Benchmarks Track Paper `ConceptFactory: Facilitate 3D Object Knowledge Annotation with Object Conceptualization`.

## Data and Template Implementations

### Dependencies

The templates are developed with the following dependencies:

- NumPy v 1.23.5
- Trimesh v 3.23.5
- Open3D v 0.17.0

### Data

The conceptualization results of objects can be found [here](https://drive.google.com/drive/folders/18fTrisH-9psUWRe8zdt4dPAyv1K4twxz).The data are organized by object categories, please download respective files (`{CATEGORY_NAME}_conceptualization.pkl`) as needed and rename the file to `conceptualization.pkl` under respective folders. 

[**Optional**] If you want to visualize both the actual object and the conceptualization results simultaneously, you can download the data from the corresponding dataset (the data IDs match one-to-one with the original dataset). Then, place the downloaded OBJ files in the `object_model` folder (you will need to create this folder yourself).

### Files

In folder `code` we provide implementations of both geometry and concept templates. The codes are organized by object categories, with each sub-folder named as `{CATEGORY_NAME}` and can be viewed as an independent module. Each of the sub-folder (after data extraction from the previous step) contains the following file/folders:

- [**Optional**] `object_model`: a folder contains the original 3D objects we have conceptualized, each object file is named as `{MODEL_ID}.obj`, where `MODEL_ID` follows the ID from the object's source, *e.g.* ShapeNet. (you will need to create this folder yourself)

- `conceptualization.pkl`: (from GoogleDrive) a list of conceptualization results.

- `base_template.py`: the definitions of the parent classes for the templates (copied across categories).

- `geometry_template.py`: the definitions of all geometry templates, which are used to construct concept templates (copied across categories).

- `concept_template.py`: the definitions of all concept templates involved in conceptualizing objects from the category.

- `utils.py`: utility functions (copied across categories).

-  `knowledge_definition.py`: the definitions of various types of knowledge organized by concepts involved in conceptualizing objects from the category. 

-  `knowledge_utils.py`: utility functions and parameters for knowledge annotation.

-  `visualize.py`: script for visualizing both the actual object and the conceptualization results (if the `object_model` folder is missing, only the conceptualization results will be displayed). The command `python visualize.py` sequentially renders object-conceptualization pairs in an Open3D window.

-  `visualize_knowledge.py`: script for visualizing the annotated knowledge on the point cloud of the objects. 

## Code Analysis & Extensions

New concept templates can be easily extended as needed. 
We have provided a detailed tutorial to achieve this in [tutorial/extension_tutorial.md](tutorial/extension_tutorial.md), which also serves
as an analysis of our current implementations.

# Citation

TBA
