# CT-TIQUA V4.0

Computed Tomography based Traumatic Brain Injury Quantification

This software takes as input a CT-scan of a head injured patient and returns: 
* 1/ the segmentation of 7 types of lesions typical from TBI
* 2/ a structural atlas dividing the input brain in 52 zones
* 3/ a vascular atlas dividing the input brain in 32 zones (from Liu, CF., Hsu, J., Xu, X. et al. Digital 3D Brain MRI Arterial Territories Atlas. Sci Data 10, 74 (2023). https://doi.org/10.1038/s41597-022-01923-0)
* 3/ Two csv files containing the volume in mm3 of all type on lesion in all zone of the brain (one for each atlas)
* 4/ A csv file containing the ventricles volume (using TotalSegmentator Wasserthal, J., Breit, H. C., Meyer, M. T., Pradella, M., Hinck, D., Sauter, A. W., ... & Segeroth, M. (2023). TotalSegmentator: robust segmentation of 104 anatomic structures in CT images. Radiology: Artificial Intelligence, 5(5).)
 
It was developed by Clément Brossard and Benjamin Lemasson (benjamin.lemasson@univ-grenoble-alpes.fr). 
Brossard C, Grèze J, de Busschère JA, Attyé A, Richard M, Tornior FD, Acquitter C, Payen JF, Barbier EL, Bouzat P, Lemasson B. 
Prediction of therapeutic intensity level from automatic multiclass segmentation of traumatic brain injury lesions on CT-scans. 
Sci Rep. 2023 Nov 17;13(1):20155. doi: 10.1038/s41598-023-46945-9. PMID: 37978266; PMCID: PMC10656472.

If you use our software, please cite our work! 
This software is not available for commercial purposes.

To use this software :
* 1/ Install docker : https://docs.docker.com/get-docker/

* 2/ Get the path of the file you want to process (we are going to call it /Path/input/data, and the file input_file.nii), and the name of the output folder (output_folder).

* 3/ Download the docker image containing all the tools needed by CT-TIQUA and run the image:
These two steps are achieved by the command:

`docker run --entrypoint=/bin/sh --rm -v /Path/input/data:/Path/input/data -w /Path/input/data -it --ipc=host nifmgin/ct-tiqua:3.2`

Details of the docker image are available here : https://hub.docker.com/repository/docker/nifmgin/ct-tiqua/general

* 4/ You are now in the interactive mode of the docker image (it can be assessed by the `#` at the beginning of your commandline). Now we can execute CT-TIQUA in the docker image:

`ct-tiqua --input input_file.nii --output output_folder`

* 5/ After the computation, you can exit the interactive mode of the docker image by typing `exit`

* 6/ You can delete all your local docker images and free the disk space by typing

`docker system prune -a -f`
