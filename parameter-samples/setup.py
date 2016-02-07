#!/usr/bin/env python

#
#   write_option_file
#
def write_option_file(option, u, file):
    # option text
    opttext = ""
    opttext += "-Metos3DDebugLevel                                  0" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DGeometryType                                Profile" + "\n"
    opttext += "-Metos3DProfileInputDirectory                       data/TMM/2.8/Geometry/" + "\n"
    opttext += "-Metos3DProfileIndexStartFile                       gStartIndices.bin" + "\n"
    opttext += "-Metos3DProfileIndexEndFile                         gEndIndices.bin" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DTracerCount                                 2" + "\n"
    opttext += "-Metos3DTracerInitValue                             2.17e+0,1.e-4" + "\n"
    opttext += "-Metos3DTracerOutputDirectory                       work/" + "\n"
    opttext += "-Metos3DTracerOutputFile                            PO4.petsc,DOP.petsc" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DParameterCount                              7" + "\n"
    opttext += "-Metos3DParameterValue                              " + ",".join(["%.16e" % ui for ui in u]) + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DBoundaryConditionCount                      2" + "\n"
    opttext += "-Metos3DBoundaryConditionInputDirectory             data/TMM/2.8/Forcing/BoundaryCondition/" + "\n"
    opttext += "-Metos3DBoundaryConditionName                       Latitude,IceCover" + "\n"
    opttext += "-Metos3DLatitudeCount                               1" + "\n"
    opttext += "-Metos3DLatitudeFileFormat                          latitude.petsc" + "\n"
    opttext += "-Metos3DIceCoverCount                               12" + "\n"
    opttext += "-Metos3DIceCoverFileFormat                          fice_$02d.petsc" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DDomainConditionCount                        2" + "\n"
    opttext += "-Metos3DDomainConditionInputDirectory               data/TMM/2.8/Forcing/DomainCondition/" + "\n"
    opttext += "-Metos3DDomainConditionName                         LayerDepth,LayerHeight" + "\n"
    opttext += "-Metos3DLayerDepthCount                             1" + "\n"
    opttext += "-Metos3DLayerDepthFileFormat                        z.petsc" + "\n"
    opttext += "-Metos3DLayerHeightCount                            1" + "\n"
    opttext += "-Metos3DLayerHeightFileFormat                       dz.petsc" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DTransportType                               Matrix" + "\n"
    opttext += "-Metos3DMatrixInputDirectory                        data/TMM/2.8/Transport/Matrix5_4/1dt/" + "\n"
    opttext += "-Metos3DMatrixCount                                 12" + "\n"
    opttext += "-Metos3DMatrixExplicitFileFormat                    Ae_$02d.petsc" + "\n"
    opttext += "-Metos3DMatrixImplicitFileFormat                    Ai_$02d.petsc" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DTimeStepStart                               0.0" + "\n"
    opttext += "-Metos3DTimeStepCount                               2880" + "\n"
    opttext += "-Metos3DTimeStep                                    0.0003472222222222" + "\n"
    opttext += "" + "\n"
    opttext += "-Metos3DSolverType                                  Newton" + "\n"
    opttext += "-Metos3DNewton_snes_type                            ls" + "\n"
    opttext += "-Metos3DNewton_snes_atol                            9.008823302130e-04" + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew" + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew_gamma                    1.0" + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew_alpha                    1.2" + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew_alpha2                   1.2" + "\n"
    opttext += "-Metos3DNewton_snes_monitor" + "\n"
    opttext += "-Metos3DNewton_ksp_type                             gmres" + "\n"
    opttext += "-Metos3DNewton_ksp_atol                             9.008823302130e-04" + "\n"
    opttext += "-Metos3DNewton_ksp_gmres_modifiedgramschmidt" + "\n"
    opttext += "-Metos3DNewton_ksp_monitor" + "\n"
    # write file
    f = open(file, "w")
    f.write(opttext)
    f.close()

#
#   write_job_file
#
def write_job_file(job, file):
    # job text
    jobtext = ""
    jobtext += "#PBS -N sample." + job + "\n"
    jobtext += "#PBS -b 4" + "\n"
    jobtext += "#PBS -l cpunum_job=16" + "\n"
    jobtext += "#PBS -l elapstim_req=01:00:00" + "\n"
    jobtext += "#PBS -l memsz_job=32gb" + "\n"
    jobtext += "#PBS -T mpich" + "\n"
    jobtext += "#PBS -j o" + "\n"
    jobtext += "#PBS -o work/metos3d.nesh-fe.sample." + job + ".out.txt" + "\n"
    jobtext += "#PBS -q clfocean" + "\n"
    jobtext += "cd $PBS_O_WORKDIR" + "\n"
    jobtext += ". job/de.uni-kiel.rz.nesh-fe.petsc-3.3-p5.opt.txt" + "\n"
    jobtext += "mpirun -r ssh $NQSII_MPIOPTS -n 64 ./metos3d-simpack-MITgcm-PO4-DOP.exe option/sample." + job + ".option.txt" + "\n"
    jobtext += "qstat -f ${PBS_JOBID/0:/}" + "\n"
    # write file
    f = open(file, "w")
    f.write(jobtext)
    f.close()

#
#   create_job_nesh_fe
#
def create_job_nesh_fe(samples):
    import numpy as np
    print "# Creating job files ..."
    # shape
    (n ,m) = samples.shape
    # generate m job files
    for i in range(0, m):
        job = "%02d" % i
        file = "job/sample." + job + ".job.txt"
#        print job, file
        write_job_file(job, file)

#
#   create_option
#
def create_option(samples):
    import numpy as np
    print "# Creating option files ..."
    # shape
    (n ,m) = samples.shape
    # generate m option files
    for i in range(0, m):
        option = "%02d" % i
        file = "option/sample." + option + ".option.txt"
#        print option, file
        write_option_file(option, samples[:,i], file)

#
#   create_job_option_nesh_fe
#
def create_job_option_nesh_fe(samples):
    create_job_nesh_fe(samples)
    create_option(samples)

#
#   read_samples
#
def read_samples():
    import numpy as np
    print "# Reading samples ..."
    # read
    fid = open('LHS-samples.bin', 'rb')
    m   = np.fromfile(fid, dtype = '>i4', count = 1)
    n   = np.fromfile(fid, dtype = '>i4', count = 1)
    x   = np.fromfile(fid, dtype = '>f8', count = m * n)
    fid.close()
    # reshape
    x = np.reshape(x, (n, m))
    return x

#
#   main
#
if __name__ == '__main__':
    samples = read_samples()
    create_job_option_nesh_fe(samples)
