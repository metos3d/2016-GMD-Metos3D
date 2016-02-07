#!/usr/bin/env python

#
#   write_job_file
#
def write_job_file(gamma, alpha):
    gstr = "{0:1.1f}".format(gamma)
    astr = "{0:.5f}".format(alpha)
    filltext = gstr + "-" + astr
    # job file
    jobfile = "job/NK-" + filltext + ".nesh-fe.job.txt"
    # job text
    jobtext = ""
    jobtext += "#PBS -N NK-" + filltext + "\n"
    jobtext += "#PBS -j o" + "\n"
    jobtext += "#PBS -o work/NK-" + filltext + ".out.txt" + "\n"
    jobtext += "#PBS -b 4" + "\n"
    jobtext += "#PBS -l cpunum_job=16" + "\n"
    jobtext += "#PBS -l elapstim_req=02:00:00" + "\n"
    jobtext += "#PBS -l memsz_job=32gb" + "\n"
    jobtext += "#PBS -T mpich" + "\n"
    jobtext += "#PBS -q clmedium" + "\n"
    jobtext += "cd $PBS_O_WORKDIR" + "\n"
    jobtext += ". job/de.uni-kiel.rz.nesh-fe.petsc-3.3-p5.opt.txt" + "\n"
    jobtext += "mpirun -n 64 -r ssh $NQSII_MPIOPTS ./metos3d-simpack-MITgcm-PO4-DOP.exe option/NK-" + filltext + ".option.txt" + "\n"
    jobtext += "qstat -f ${PBS_JOBID/0:/}" + "\n"
    # store
    f = open(jobfile, "w")
    f.write(jobtext)
    f.close()

#
#   write_option_file
#
def write_option_file(gamma, alpha):
    gstr = "{0:1.1f}".format(gamma)
    astr = "{0:.5f}".format(alpha)
    filltext = gstr + "-" + astr
    # opt file
    optfile = "option/NK-" + filltext + ".option.txt"
    # opt text
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
    opttext += "-Metos3DParameterValue                              0.5,2.0,0.67,0.5,30.0,0.02,0.858" + "\n"
    # boundary data
    opttext += "" + "\n"
    opttext += "-Metos3DBoundaryConditionCount                      2" + "\n"
    opttext += "-Metos3DBoundaryConditionInputDirectory             data/TMM/2.8/Forcing/BoundaryCondition/" + "\n"
    opttext += "-Metos3DBoundaryConditionName                       Latitude,IceCover" + "\n"
    opttext += "-Metos3DLatitudeCount                               1" + "\n"
    opttext += "-Metos3DLatitudeFileFormat                          latitude.petsc" + "\n"
    opttext += "-Metos3DIceCoverCount                               12" + "\n"
    opttext += "-Metos3DIceCoverFileFormat                          fice_$02d.petsc" + "\n"
    # domain data
    opttext += "" + "\n"
    opttext += "-Metos3DDomainConditionCount                        2" + "\n"
    opttext += "-Metos3DDomainConditionInputDirectory               data/TMM/2.8/Forcing/DomainCondition/" + "\n"
    opttext += "-Metos3DDomainConditionName                         LayerDepth,LayerHeight" + "\n"
    opttext += "-Metos3DLayerDepthCount                             1" + "\n"
    opttext += "-Metos3DLayerDepthFileFormat                        z.petsc" + "\n"
    opttext += "-Metos3DLayerHeightCount                            1" + "\n"
    opttext += "-Metos3DLayerHeightFileFormat                       dz.petsc" + "\n"
    # transport
    opttext += "" + "\n"
    opttext += "-Metos3DTransportType                               Matrix" + "\n"
    opttext += "-Metos3DMatrixInputDirectory                        data/TMM/2.8/Transport/Matrix5_4/1dt/" + "\n"
    opttext += "-Metos3DMatrixCount                                 12" + "\n"
    opttext += "-Metos3DMatrixExplicitFileFormat                    Ae_$02d.petsc" + "\n"
    opttext += "-Metos3DMatrixImplicitFileFormat                    Ai_$02d.petsc" + "\n"
    # time step
    opttext += "" + "\n"
    opttext += "-Metos3DTimeStepStart                               0.0" + "\n"
    opttext += "-Metos3DTimeStepCount                               2880" + "\n"
    opttext += "-Metos3DTimeStep                                    0.0003472222222222" + "\n"
    # solver
    opttext += "-Metos3DSolverType                                  Newton" + "\n"
    opttext += "-Metos3DNewton_snes_type                            ls" + "\n"
    opttext += "-Metos3DNewton_snes_atol                            9.008823302130e-04" + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew" + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew_gamma                    " + gstr + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew_alpha                    " + astr + "\n"
    opttext += "-Metos3DNewton_snes_ksp_ew_alpha2                   " + astr + "\n"
    opttext += "-Metos3DNewton_snes_monitor" + "\n"
    opttext += "-Metos3DNewton_ksp_type                             gmres" + "\n"
    opttext += "-Metos3DNewton_ksp_atol                             9.008823302130e-04" + "\n"
    opttext += "-Metos3DNewton_ksp_max_it                           300" + "\n"
    opttext += "-Metos3DNewton_ksp_gmres_modifiedgramschmidt" + "\n"
    opttext += "-Metos3DNewton_ksp_monitor" + "\n"
    # store
    f = open(optfile, "w")
    f.write(opttext)
    f.close()

#
#   create_job_option
#
def create_job_option():
    import math as m
    # ew params
    gammas = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5]
    alphas = [(1 + m.sqrt(5))/2, 1.5, 1.4, 1.3, 1.2, 1.1]
    # files
    for gamma in gammas:
        for alpha in alphas:
            write_job_file(gamma, alpha)
            write_option_file(gamma, alpha)

#
#   main
#
if __name__ == "__main__":
    create_job_option()

