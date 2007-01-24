import sys,os,string

if __name__=="__main__":
    if( len(sys.argv) >= 2 ):
        paths = sys.argv[1:]
    for path in paths:
        outpath = os.path.splitext( path )[0]+".py"
        infile = open(path,"r")
        outfile = open(outpath,"w")
        outfile.write('"""\n')
        outfile.write(outpath+"\n")
        outfile.write("Automatically generated.\n")
        outfile.write('"""\n')
        lines = infile.readlines()
        for line in lines:
            if ( line.find("::") > 0 ):
                line = line.split("::")[1]
                line = line.split("!")[0]
                pieces = line.split("=")
                outfile.write( "%s = %s\n" %(pieces[0],pieces[1]) )
            
