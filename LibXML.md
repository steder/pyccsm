# Obtaining: #

[ftp://xmlsoft.org/libxml2/](ftp://xmlsoft.org/libxml2/)
```
lslogin1% wget ftp://xmlsoft.org/libxml2/libxml2-2.6.11.tar.gz
lslogin1% pwd
/san/hpc/A-ig2/steder/Sandbox
lslogin1% ls
babel-1.0.2         chasm                  MCT               Numeric-24.2.installed  numpy-1.0.1.tar.gz
babel-1.0.2.tar.gz  chasm_1.4.RC1.tar.gz   MCTnodata.tar.gz  Numeric-24.2.tar.gz     Python-2.5
babel.tar.gz        libxml2-2.6.11.tar.gz  Numeric-24.2      numpy-1.0.1             Python-2.5.tgz
lslogin1% tar -zxf libxml2-2.6.11.tar.gz 
lslogin1% cd libxml2-2.6.11
lslogin1% ls
acconfig.h                    install-sh                    testThreads.c
acinclude.m4                  legacy.c                      testThreadsWin32.c
aclocal.m4                    libxml-2.0.pc.in              testURI.c
AUTHORS                       libxml-2.0-uninstalled.pc.in  testXPath.c
c14n.c                        libxml2.spec                  threads.c
catalog.c                     libxml.3                      TODO
ChangeLog                     libxml.h                      TODO_SCHEMAS
check-relaxng-test-suite2.py  libxml.m4                     tree.c
check-relaxng-test-suite.py   libxml.spec.in                trio.c
check-xinclude-test-suite.py  list.c                        triodef.h
check-xml-test-suite.py       ltmain.sh                     trio.h
check-xsddata-test-suite.py   macos                         trionan.c
chvalid.c                     Makefile.am                   trionan.h
config.guess                  Makefile.in                   triop.h
config.h.in                   missing                       triostr.c
config.sub                    mkinstalldirs                 triostr.h
configure                     nanoftp.c                     uri.c
configure.in                  nanohttp.c                    valid.c
COPYING                       NEWS                          vms
Copyright                     parser.c                      win32
dbgenattr.pl                  parserInternals.c             xinclude.c
dbgen.pl                      pattern.c                     xlink.c
debugXML.c                    python                        xml2-config.1
depcomp                       README                        xml2-config.in
dict.c                        regressions.py                xml2Conf.sh.in
doc                           regressions.xml               xmlcatalog.c
DOCBparser.c                  relaxng.c                     xmlIO.c
elfgcchack.h                  result                        xmllint.c
encoding.c                    SAX2.c                        xmlmemory.c
entities.c                    SAX.c                         xmlreader.c
error.c                       test                          xmlregexp.c
example                       testAutomata.c                xmlsave.c
genUnicode.py                 testC14N.c                    xmlschemas.c
globals.c                     testHTML.c                    xmlschemastypes.c
hash.c                        testReader.c                  xmlstring.c
HTMLparser.c                  testRegexp.c                  xmlunicode.c
HTMLtree.c                    testRelax.c                   xmlwriter.c
include                       testSAX.c                     xpath.c
INSTALL                       testSchemas.c                 xpointer.c
lslogin1% 
```

# Config #

```
$ ./configure --prefix=$PYMCT_ROOT --with-python=$PYTHON
```

# Build #

```
$ make
```

# Install #
```
$ make install
```