import YAML

function get_params(fname)

    data = YAML.load(open(fname))

    key = "grid"
    xmin = data[key]["xmin"]
    xmax = data[key]["xmax"]
    ymin = data[key]["xmin"]
    ymax = data[key]["xmax"]
    nx = data[key]["nx"]
    ny = data[key]["ny"]
    ghostwidth = data[key]["ghostwidth"]

    key = "numerics"
    itmax = data[key]["itmax"]
    cfl = data[key]["cfl"]

    key = "output"
    vtkfreq = data[key]["vtkfreq"]
    
    return Int64[nx, ny, itmax, vtkfreq, ghostwidth], Float64[xmin, xmax, ymin, ymax, cfl]
end
