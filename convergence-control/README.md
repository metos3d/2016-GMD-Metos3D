#### Setup:

    $> . job/de.uni-kiel.rz.nesh-fe.petsc-3.3-p5.opt.txt
    $> metos3d simpack MITgcm-PO4-DOP

#### Run:

    $> for i in $(ls job/NK-*); do qsub $i; done;
