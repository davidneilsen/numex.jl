import YAML
include("pars.jl")

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
    dtype = data[key]["dtype"]
    dissipation = data[key]["dissipation"]

    key = "output"
    vtkfreq = data[key]["vtkfreq"]
    asciifreq = data[key]["asciifreq"]

    key = "initdata"
    idtype = data[key]["idtype"]
    idwg_bzamp = data[key]["idwg_bzamp"]
    idwg_omega = data[key]["idwg_omega"]
    idwg_m_mode = data[key]["idwg_m_mode"]
    idwg_n_mode = data[key]["idwg_n_mode"]

    if dtype == 1
        dtype = min(2, dtype)
    elseif dtype == 2
        dtype = min(3, dtype)
    end

    pint = zeros(Int64, n_params_int)
    pint[p_nx] = nx
    pint[p_ny] = ny
    pint[p_itmax] = itmax
    pint[p_vtkfreq] = vtkfreq
    pint[p_ghostwidth] = ghostwidth
    pint[p_dtype] = dtype
    pint[p_idtype] = idtype
    pint[p_idwg_m_mode] = idwg_m_mode
    pint[p_idwg_n_mode] = idwg_n_mode
    pint[p_asciifreq] = asciifreq
    pint[p_dissipation] = dissipation
    
    preal = zeros(Float64, n_params_double)
    preal[p_xmin] = xmin
    preal[p_xmax] = xmax
    preal[p_ymin] = ymin
    preal[p_ymax] = ymax
    preal[p_cfl] = cfl
    preal[p_idwg_bzamp] = idwg_bzamp

    return pint, preal
end

function buildpars(pint, preal)

    pars = Dict("itmax"=>pint[p_itmax], "vtkfreq"=>pint[p_vtkfreq], "asciifreq"=> pint[p_asciifreq], "ghostwidth"=>pint[p_ghostwidth], "dtype"=>pint[p_dtype], "idtype"=>pint[p_idtype], "idwg_m_mode"=>pint[p_idwg_m_mode], "idwg_n_mode"=>pint[p_idwg_n_mode], "cfl"=>preal[p_cfl], "idwg_bzamp"=>preal[p_idwg_bzamp], "dissipation"=>pint[p_dissipation])

    return pars
end
