using MPI

MPI.Init()
comm = MPI.COMM_WORLD
nnodes = MPI.Comm_size(comm)
ndims = 3
dims = MPI.Dims_create(nnodes, [0,0,0])
MPI.Finalize()

