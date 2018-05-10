from py4j.java_gateway import JavaGateway

"""
    cal = Filename of calibration signal
    out = Filename of sensor output to calibration
    resp = Filename of response to load in
    embd = Set this to true if response is embedded in the jar
    start = Lower bound of time region of interest, in milliseconds since epoch
    end = Upper bound of region of interest, in milliseconds since epoch
    lf = True if the data under examination is a low-frequency cal
"""
def getCalc(cal, out, resp, embd, start, end, lf):
    gateway = JavaGateway()
    exp = 
	gateway.entry_point.populateDataAndRun(cal, out, resp, embd, start, end, lf)
    zerosIn = exp.getInitialZeros();
    zerosOut = exp.getFitZeros();
    polesIn = exp.getInitialPoles();
    polesOut = exp.getFitPoles();
    initResid = exp.getInitResidual();
    fitResid = exp.getFitResidual();
    freqs = exp.getFreqList();
    magCurves = exp.getAmplitudesAsArrays();
    phaseCurves = exp.getPhasesAsArrays();

