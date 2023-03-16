function nint(x)
    if (x >= 0.0)
        n = round(Int, x)
    elseif
        n = round(Int, x)
    end
    return n
end

function cluster!(gh, dims, rboxcount)

    boxcount = 0
    shp0 = gh.shp0

    minx = 1
    maxx = shp0[1]
    miny = 1
    maxy = shp0[2]

    lengthi = nint((maxx - minx)/dims[1])
    lengthj = nint((maxy - miny)/dims[2])

    numx = dims[1]
    numy = dims[2]

    if ( (numx > 1 || numy > 1) && (gh.maxnum >= (boxcount + numx*numy)) )
        for j = 1:numy
            for i = 1:numx
                boxcount += 1

                gh.ibox[boxcount,1,1] = minx + (i-1)*lengthi
                gh.ibox[boxcount,1,2] = maxx + (i)*lengthi
                gh.ibox[boxcount,2,1] = miny + (j-1)*lengthj
                gh.ibox[boxcount,2,2] = maxy + (j)*lengthj

                gh.gibox[boxcount,1,1] = minx + (i-1)*lengthi
                gh.gibox[boxcount,1,2] = minx + (i)*lengthi
                gh.gibox[boxcount,2,1] = miny + (i-1)*lengthi
                gh.gibox[boxcount,2,2] = miny + (i)*lengthi

                if i == 1 gh.ibox[boxcount,1,1] = minx end
                if j == 1 gh.ibox[boxcount,2,1] = miny end
                if i == numx gh.ibox[boxcount,1,2] = maxx end
                if j == numy gh.ibox[boxcount,2,2] = maxy end

                if i == 1 gh.gibox[boxcount,1,1] = minx end
                if j == 1 gh.gibox[boxcount,2,1] = miny end
                if i == numx gh.gibox[boxcount,1,2] = maxx end
                if j == numy gh.gibox[boxcount,2,2] = maxy end

                # Add ghostzones
                if gh.ibox[boxcount,1,1] != minx gh.ibox[boxcount,1,1] -= gh.ghostwidth end
                if gh.ibox[boxcount,1,2] != maxx gh.ibox[boxcount,1,1] += gh.ghostwidth end
                if gh.ibox[boxcount,2,1] != miny gh.ibox[boxcount,1,1] -= gh.ghostwidth end
                if gh.ibox[boxcount,2,2] != maxy gh.ibox[boxcount,1,1] += gh.ghostwidth end

            end
        end
    else
        boxcount += 1

        gh.ibox[boxcount,1,1] = minx
        gh.ibox[boxcount,1,2] = maxx
        gh.ibox[boxcount,2,1] = miny
        gh.ibox[boxcount,2,2] = maxy
        
        gh.gibox[boxcount,1,1] = minx
        gh.gibox[boxcount,1,2] = maxx
        gh.gibox[boxcount,2,1] = miny
        gh.gibox[boxcount,2,2] = maxy
    end

    rboxcount = boxcount
end




function grid_sync( nfields)

    return 0 if gh.numRanks == 1
    
    send_buffer = Array{Array{Float64,1},1}(undef, gh.comm_count)
    recv_buffer = Array{Array{Float64,1},1}(undef, gh.comm_count)

    for i = 1:gh.comm_count
        send_buffer[i] = Array{Float64,1}(undef, numfields*gh.length1[i])
        recv_buffer[i] = Array{Float64,1}(undef, numfields*gh.length2[i])
    end

    
end