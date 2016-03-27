from gravpy import nfwf

print "running nfwf module"
print nfwf.single_eval(0.0,0.0,(1,0,0,0.0001,0,0.1))[1:]
print nfwf.single_eval(0.1,0.0,(1,0,0,0.0001,0,0.1))[1:]

