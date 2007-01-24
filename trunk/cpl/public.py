import os,sys

if __name__=="__main__":
    if ( len( sys.argv ) >= 2 ):
        lines = open( sys.argv[1], "r" ).readlines()
        publics = []
        for line in lines:
            if (line.find("public")>=0):
                publics.append(line)
        outfile = open( sys.argv[1] + ".public","w")
        for line in publics:
            outfile.write(line)
    else:
        print "Usage:"
        print "\t./public.py (filename)"
