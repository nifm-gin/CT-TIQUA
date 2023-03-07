# CT-TIQUA

Computed Tomography based Traumatic Brain Injury Quantification

This software takes as input a CT-scan of a head injured patient and returns: 
* 1/ the segmentation of 7 types of lesions typical from TBI
* 2/ a structural atlas dividing the input brain in 10 zones
* 3/ a csv file containing the volume in mm3 of all type on lesion in all zone of the brain 
 
It was developed by Clément Brossard and Benjamin Lemasson (benjamin.lemasson@univ-grenoble-alpes.fr). We are currently writing a scientific article descibing the whole process. If you use our software, please cite our work! Part of this repository is taken and modified from the repository https://github.com/biomedia-mira/blast-ct , so please consider citing their work too if relevant.

To use this software :
* 1/ Install docker : https://docs.docker.com/get-docker/

* 2/ Get the path of the file you want to process (we are going to call it /Path/input/data, and the file input_file.nii), and the name of the output folder (output_folder).

* 3/ Download the docker image containing all the tools needed by CT-TIQUA and run the image:
These two steps are achieved by the command:

`docker run --entrypoint=/bin/sh --rm -v /Path/input/data:/Path/input/data -w /Path/input/data -it --ipc=host brosscle/ct-tiqua:1.3`

Details of the docker image are available here : https://hub.docker.com/repository/docker/brosscle/ct-tiqua/general

* 4/ You are now in the interactive mode of the docker image (it can be assessed by the `#` at the beginning of your commandline). Now we can execute CT-TIQUA in the docker image:

`ct-tiqua --input input_file.nii --output output_folder --ensemble`

* 5/ After the computation, you can exit the interactive mode of the docker image by typing `exit`

* 6/ You can delete all your local docker images and free the disk space by typing

`docker system prune -a -f`
