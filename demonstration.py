#!/opt/local/bin/python2.7

import BBmethodsNEURON as bbm

data = bbm.runAC_vLS(AC = 0.4, vLS = 7.0, Istim=12.0, LoadFromPrevious=False)
bbm.plot(data, xmin=None, xmax=None, InteractiveFigure=False)


