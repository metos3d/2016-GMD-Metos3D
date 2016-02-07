### Latin Hypercube Sampling

Generated samples using MATLAB's `lhsdesign` routine:

    m  = 100;
    n  = 7;
    xx = lhsdesign(m, n);

    lb = [0.25,   1.5, 0.05, 0.25, 10.0, 0.01, 0.7];
    ub = [0.75, 200.0, 0.95, 1.50, 50.0, 0.05, 1.5];

    x  = bsxfun(@plus, lb, bsxfun(@times, xx, (ub-lb)))

#### Written to disk using MATLAB:

    fid = fopen('LHS-samples.bin', 'wb', 'ieee-be')
    fwrite(fid, m, 'integer*4')
    fwrite(fid, n, 'integer*4')
    fwrite(fid, x, 'real*8')
    fclose(fid)

#### Read in using Python:

    import numpy as np

    fid = open('LHS-samples.bin', 'rb')
    m   = np.fromfile(fid, dtype = '>i4', count = 1)
    n   = np.fromfile(fid, dtype = '>i4', count = 1)
    x   = np.fromfile(fid, dtype = '>f8', count = m * n)
    fid.close()

#### Reshaped:

    x = np.reshape(x, (n, m))

### Setup:

    $> . job/de.uni-kiel.rz.nesh-fe.petsc-3.3-p5.opt.txt
    $> metos3d simpack MITgcm-PO4-DOP

### Run:

    $> for i in $(ls job/sample.*); do qsub $i; done;
