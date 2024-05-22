#=-------------------------------------------------------------------
 =
 =  This routine writes the data to an ASCII curve file that can be 
 =  read by VisIT.  The curve format is only used by VisIT, I think.
 =
 =  fname = the name of the file without the .curve extension
 =  names = the function names in an array. (array of strings)
 =
 =-------------------------------------------------------------------=#


function output_curve(fname, names, time, fields)
    nx = fields.grid.nx
    x = fields.grid.x
    neqs = fields.neqs
    u = fields.u

    filename = @sprintf("%s.curve", fname)
    fp = open(filename, "w")
        
    @printf(fp,"#TIME %14.7e\n", time)
    for m = 1:neqs
        @printf(fp, "# %s\n",names[m]);
        for i = 1:nx
            @printf(fp,"%14.7e    %14.7e\n",x[i],u[m][i])
        end
    end
    close(fp)
end
    
