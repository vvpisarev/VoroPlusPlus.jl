function VoronoiCell(xlo, xhi, ylo, yhi, zlo, zhi)
    xmin, xmax, ymin, ymax, zmin, zmax = Float64.((xlo, xhi, ylo, yhi, zlo, zhi))
    v = VoronoiCell()
    init!(v, xmin, xmax, ymin, ymax, zmin, zmax)
    return v
end
