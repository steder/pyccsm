import traceback

modules = ['bundle','comm','const','control','contract','domain','error',
           'fields','infobuffer']
def test():
    for module in modules:
        try:
            exec( "import %s" % module )
        except:
            traceback.print_exc()

if __name__=="__main__":
    test()
