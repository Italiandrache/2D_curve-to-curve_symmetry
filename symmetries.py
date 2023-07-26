import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from itertools import islice
import multiprocessing
import time


def notReal(x_t, tPy, t, tNum):
    t = sp.Symbol('t', real = True) if tPy == "t" else sp.Symbol('q', real = True) 
    if t in sp.sympify(x_t).free_symbols:
        if not x_t.subs(t, tNum).is_real:
            return True
    return False

def returnValue(x_t, tPy, t, tNum):
    t = sp.Symbol('t', real = True) if tPy == "t" else sp.Symbol('q', real = True)
    if t in sp.sympify(x_t).free_symbols:
        return x_t.subs(t, tNum)
    return x_t

def returnValueNoAss(x_t, t, tNum):
    if t in sp.sympify(x_t).free_symbols:
        return x_t.subs(t, tNum)
    return x_t

def curve(x_t, y_t, tPy, t, tNum):
    if notReal(x_t, tPy, t, tNum) or notReal(y_t, tPy, t, tNum):
        return 0, 0, False
    return returnValue(x_t, tPy, t, tNum), returnValue(y_t, tPy, t, tNum), True

def getCoeffPerp(x_t, y_t, xPrime_t1, yPrime_t1, tPy, t, t1, tNum):
    yPrime = returnValueNoAss(yPrime_t1, t1, tNum)
    xPrime = returnValueNoAss(xPrime_t1, t1, tNum)
    if not (yPrime == sp.sympify("nan") or yPrime == sp.sympify("+oo") or yPrime == sp.sympify("-oo") or yPrime == sp.sympify("zoo")): #derivative at tNum is a real number. tNum is always real so no chanche of having a derivative with non-zero immaginary part
        yPrime_t = sp.diff(y_t, t) if isinstance(y_t, sp.Basic) else 0 #calculate the derivative with the real=True assumption, such that it returns just a number
        yPrime = returnValue(yPrime_t, tPy, t, tNum)
    if not (xPrime == sp.sympify("nan") or xPrime == sp.sympify("+oo") or xPrime == sp.sympify("-oo") or xPrime == sp.sympify("zoo")):
        xPrime_t = sp.diff(x_t, t) if isinstance(x_t, sp.Basic) else 0
        xPrime = returnValue(xPrime_t, tPy, t, tNum)
    
    if yPrime == sp.sympify("nan") or xPrime == sp.sympify("nan") or (not yPrime.is_real and not xPrime.is_real) or (yPrime == 0 and xPrime == 0):
        return np.nan
    elif yPrime != 0 and xPrime != 0:
        return -xPrime/yPrime #-1/coeff tan
    elif xPrime == 0:
        return 0
    else: #yPrime == 0
        return np.inf

def segment(xToBeMirrored_t, yToBeMirrored_t, tPy, t, tRange, tIntervals = ()):
    xSegmentList = []
    ySegmentList = []
    tNum_1List = []
    tNum_2List = []
    q = sp.Symbol('q', real = True)
    for i in range(0, len(tRange)-1):
        tNum_1 = tRange[i]
        tNum_2 = tRange[i+1]
        xToBeMirrored_1, yToBeMirrored_1, real_1 = curve(xToBeMirrored_t, yToBeMirrored_t, tPy, t, tNum_1)
        xToBeMirrored_2, yToBeMirrored_2, real_2 = curve(xToBeMirrored_t, yToBeMirrored_t, tPy, t, tNum_2)
        if not real_1 or not real_2:
            continue
        piecewiseExcep = False
        if len(tIntervals) != 0:
            tIntervalsList = []
            if isinstance(xToBeMirrored_t, sp.Piecewise) or isinstance(yToBeMirrored_t, sp.Piecewise):
                for tInterval in tIntervals:
                    tValMin, tValMax = tInterval.args[0], tInterval.args[1]
                    tValMinOpen, tValMaxOpen = tInterval.args[2], tInterval.args[3]
                    tIntervalsList.append(((tValMin, tValMinOpen), (tValMax, tValMaxOpen)))
            for k in range(len(tIntervalsList)-1):
                if sp.sympify(tNum_1) == sp.sympify(tIntervalsList[k][1][0]) and sp.sympify(tNum_1) == sp.sympify(tIntervalsList[k+1][0][0]) and (tIntervalsList[k][1][1] == False or tIntervalsList[k+1][0][1] == False):
                    if xToBeMirrored_1 != xToBeMirrored_2 or yToBeMirrored_1 != yToBeMirrored_2:
                        piecewiseExcep = True
                        break
        if piecewiseExcep:
            continue
        xSegmentList += [xToBeMirrored_1+q*(xToBeMirrored_2-xToBeMirrored_1)]
        ySegmentList += [yToBeMirrored_1+q*(yToBeMirrored_2-yToBeMirrored_1)]
        tNum_1List += [tNum_1]
        tNum_2List += [tNum_2]
    return xSegmentList, ySegmentList, tNum_1List, tNum_2List

def intersect(xMirror, yMirror, coeff, xSegmentList, ySegmentList, qNum_1List, qNum_2List, tNum, xMirror_t, yMirror_t, xToBeMirrored_q, yToBeMirrored_q, tPy, qPy):
    intersections = []
    q = sp.Symbol('q', real = True)
    t = sp.Symbol('t', real = True)
    r = sp.Symbol('r', real = True)
    for i in range(0, len(xSegmentList)-1):
        xLine_r = xMirror+r if coeff != np.inf else xMirror
        yLine_r = yMirror+r*coeff if coeff != np.inf else r

        intersection = None
        try: #could be done more efficiently with some if/else
            intersection = sp.linsolve([xSegmentList[i]-xLine_r, ySegmentList[i]-yLine_r], q, r)
            qNum, rNum = list(intersection)[0] #if the same, qNum is symbolical. If without intersections, length is 0
            if r in qNum.free_symbols or q in qNum.free_symbols:
                coin, toBeMirrored_1Tuple, toBeMirrored_2Tuple = coincident(xMirror_t, yMirror_t, t, tNum, xToBeMirrored_q, yToBeMirrored_q, q, qNum_1List[i], qNum_2List[i], coeff, tPy, qPy)
                if coin:
                    intersections += [toBeMirrored_1Tuple]
                    intersections += [toBeMirrored_2Tuple]
            elif 0 <= qNum <= 1:
                intersections += [(xSegmentList[i].subs(q, qNum), ySegmentList[i].subs(q, qNum))]
        except ValueError:
            continue
        except IndexError:
            #probably just a continue would do
            if len(list(intersection)) != 0:
                coin, toBeMirrored_1Tuple, toBeMirrored_2Tuple = coincident(xMirror_t, yMirror_t, t, tNum, xToBeMirrored_q, yToBeMirrored_q, q, qNum_1List[i], qNum_2List[i], coeff, tPy, qPy)
                if coin:
                    intersections += [toBeMirrored_1Tuple]
                    intersections += [toBeMirrored_2Tuple]
        except Exception:
            continue
    return intersections

def linIndip(xMirror, yMirror, xToBeMirrored_1, yToBeMirrored_1, xToBeMirrored_2, yToBeMirrored_2):
    a = sp.Symbol('a', real = True)
    b = sp.Symbol('b', real = True)
    dependency = list(sp.linsolve([a*(xToBeMirrored_1-xMirror)+b*(xToBeMirrored_2-xMirror), a*(yToBeMirrored_1-yMirror)+b*(yToBeMirrored_2-yMirror)], a, b))
    aNum, bNum = dependency[0]
    if isinstance(aNum, int) and isinstance(bNum, int): #maybe useless
        if aNum == 0 and bNum == 0:
            return True
    return False

def coincident(xMirror_t, yMirror_t, t, tNum, xToBeMirrored_q, yToBeMirrored_q, q, qNum_1, qNum_2, coeff, tPy, qPy):
    xMirror = returnValue(xMirror_t, tPy, t, tNum)
    yMirror = returnValue(yMirror_t, tPy, t, tNum)
    xToBeMirrored_1 = returnValue(xToBeMirrored_q, qPy, q, qNum_1)
    xToBeMirrored_2 = returnValue(xToBeMirrored_q, qPy, q, qNum_2)
    yToBeMirrored_1 = returnValue(yToBeMirrored_q, qPy, q, qNum_1)
    yToBeMirrored_2 = returnValue(yToBeMirrored_q, qPy, q, qNum_2)
    if xToBeMirrored_1 == xToBeMirrored_2:
        if not linIndip(xMirror, yMirror, xToBeMirrored_1, yToBeMirrored_1, xToBeMirrored_2, yToBeMirrored_2):
            return (True, (xToBeMirrored_1, yToBeMirrored_1), (xToBeMirrored_2, yToBeMirrored_2))
    elif abs(coeff - (yToBeMirrored_2 - yToBeMirrored_1)/(xToBeMirrored_2 - xToBeMirrored_1)) <= 0.00001 and not linIndip(xMirror, yMirror, xToBeMirrored_1, yToBeMirrored_1, xToBeMirrored_2, yToBeMirrored_2):
        return (True, (xToBeMirrored_1, yToBeMirrored_1), (xToBeMirrored_2, yToBeMirrored_2))
    return (False, (0, 0), (0, 0))

def mirror(xSegmentList, ySegmentList, qNum_1List, qNum_2List, xMirror_t, yMirror_t, xToBeMirrored_q, yToBeMirrored_q, t, tRange, tPy, qPy, currentProcess, nProcesses, mirroredShared):
    tRange = islice(tRange, currentProcess, len(tRange), nProcesses)
    t1 = sp.Symbol('t')
    xMirror_t1, yMirror_t1 = None, None
    if t in sp.sympify(xMirror_t).free_symbols:
        xMirror_t1 = xMirror_t.subs(t, t1)
    else:
        xMirror_t1 = xMirror_t
    if t in sp.sympify(yMirror_t).free_symbols:
        yMirror_t1 = yMirror_t.subs(t, t1)
    else:
        yMirror_t1 = yMirror_t
    xPrime_t1 = sp.diff(xMirror_t1, t1) if isinstance(xMirror_t1, sp.Basic) else 0 #evaluating complex derivatives to avoid sympy inaccuracies due to assumptions (eg. sign(0) = 0 instead of nan)
    yPrime_t1 = sp.diff(yMirror_t1, t1) if isinstance(yMirror_t1, sp.Basic) else 0
    for tNum in tRange:
        print(tNum)
        xMirror, yMirror, real = curve(xMirror_t, yMirror_t, tPy, t, tNum)
        if not real:
            continue
        coeff = getCoeffPerp(xMirror_t, yMirror_t, xPrime_t1, yPrime_t1, tPy, t, t1, tNum)
        match coeff:
            case np.nan:
                continue
            case _:
                pass
        intersections = intersect(xMirror, yMirror, coeff, xSegmentList, ySegmentList, qNum_1List, qNum_2List, tNum, xMirror_t, yMirror_t, xToBeMirrored_q, yToBeMirrored_q, tPy, qPy)
        for i in intersections:
            mirroredShared.append(calcSymm(xMirror, yMirror, i[0], i[1]))

def calcSymm(xMirror, yMirror, xToBeMirroredIntersection, yToBeMirroredIntersection):
    return (2*xMirror-xToBeMirroredIntersection, 2*yMirror-yToBeMirroredIntersection)

def points(x_t, y_t, tPy, t, tRange):
    x_tList = []
    y_tList = []
    for tNum in tRange:
        if notReal(x_t, tPy, t, tNum) or notReal(y_t, tPy, t, tNum):
            continue
        x_tList += [returnValue(x_t, tPy, t, tNum)]
        y_tList += [returnValue(y_t, tPy, t, tNum)]
    return x_tList, y_tList

def getAbsCurvature(x_t, y_t, tPy, t, tNum):
    t1 = sp.Symbol('t') if tPy == "t" else sp.Symbol('q')
    x_t1, y_t1 = None, None
    if t in sp.sympify(x_t).free_symbols:
        x_t1 = x_t.subs(t, t1)
    else:
        x_t1 = x_t
    if t in sp.sympify(y_t).free_symbols:
        y_t1 = y_t.subs(t, t1)
    else:
        y_t1 = y_t
    curv_t1 = sp.sympify(sp.Abs(sp.diff(x_t1, t1)*sp.diff(y_t1, t1, 2)-sp.diff(y_t1, t1)*sp.diff(x_t1, t1, 2))/(sp.diff(x_t1, t1)**2+sp.diff(y_t1, t1)**2)**(3/2)) #absolute curvature equation
    curv = returnValueNoAss(curv_t1, t1, tNum)
    if not(curv == sp.sympify("nan") or curv == sp.sympify("+oo") or curv == sp.sympify("zoo")):
        curv_t = sp.sympify(sp.Abs(sp.diff(x_t, t)*sp.diff(y_t, t, 2)-sp.diff(y_t, t)*sp.diff(x_t, t, 2))/(sp.diff(x_t, t)**2+sp.diff(y_t, t)**2)**(3/2))
        return returnValue(curv_t, tPy, t, tNum)
    return curv

def sortMixedList(mixedList, reverseBool):
    notNum_part = []
    for i in mixedList:
        if not i.is_real or i == sp.sympify("+oo") or i == sp.sympify("zoo") or i == sp.sympify("nan"):
            notNum_part.append(i)
    num_part = sorted([i for i in mixedList if i.is_real and not (i == sp.sympify("+oo") or i == sp.sympify("zoo") or i == sp.sympify("nan"))], reverse=reverseBool)
    notNum_partLen = len(notNum_part)
    if notNum_partLen != 0:
        return num_part + notNum_part, notNum_partLen
    return num_part, notNum_partLen

def getNum(curv, curvMax, numMax, numMin):
    if curv >= curvMax:
        return numMax
    num = curv/curvMax*numMax
    if num <= numMin:
        return numMin
    return num

def generateRange(rangeValuesList, variableDensities=False, x_t=None, y_t=None, tPy=None, t=None):
    for i in range(len(rangeValuesList)):
        if len(rangeValuesList[i]) == 3:
            rangeValuesList[i] = rangeValuesList[i] + (rangeValuesList[i][2],) #adding a numMin witch is equal to numMax inside the tuple, so to be able not to define numMin in tRangeValueList when we only want a fixed density for that interval
    if not variableDensities: #density constant everywhere
        tRange = np.hstack([np.linspace(rangeValues[0], rangeValues[1], num=rangeValues[2]) for rangeValues in rangeValuesList])
        for i in reversed(range(1, len(tRange))):
            if tRange[i] == tRange[i-1]:
                tRange = np.delete(tRange, i)
        return tRange
    tRange = np.array([])
    for rangeValues in rangeValuesList:
        if rangeValues[2] == rangeValues[3]: #density set to be constant in this interval. No need for below computations
            tRange = np.concatenate((tRange, np.linspace(rangeValues[0], rangeValues[1], num=rangeValues[2])))
            continue
        Delta_t, tNum = 0, rangeValues[0]
        tRangeCurv = np.linspace(rangeValues[0], rangeValues[1], num=rangeValues[2]) #crete a range in which evaluating curvature using highest allowed density
        curvMax, curvList, edgeCase = 0, [], False
        for tNumCurv in tRangeCurv:
            curvList.append(getAbsCurvature(x_t, y_t, tPy, t, tNumCurv))
        
        curvListSorted, notNumLen = sortMixedList(curvList, True) #reverse=True -> sort curve from max value to min value
        if curvListSorted[0].is_real and not curvListSorted[0] == sp.sympify("+oo") and not curvListSorted[0] == sp.sympify("zoo") and not curvListSorted[0] == sp.sympify("nan"):
            curvMax = curvListSorted[0]
        else:
            edgeCase = True #complex or infinite or undefined (NaN) curvature everywhere
        
        if edgeCase or curvListSorted[0] == curvListSorted[-1]: #density happens to be constant in this interval
            tRange = np.concatenate((tRange, tRangeCurv))
            continue

        if notNumLen != 0:
            replacedValues = 0
            for i in range(len(curvList)): #replace non real non finite or undefined curv values with curvMax
                if not curvList[i].is_real or curvList[i] == sp.sympify("+oo") or curvList[i] == sp.sympify("zoo") or curvList[i] == sp.sympify("nan"):
                    curvList[i] = curvMax
                    replacedValues += 1
                if replacedValues == notNumLen:
                    break
        
        curv = 0
        tRange = np.append(tRange, rangeValues[0])
        while tNum < rangeValues[1]:
            print(tRange.size)
            if tNum == rangeValues[0]:
                curv = curvList[0]
            else:
                curv = getAbsCurvature(x_t, y_t, tPy, t, tNum)
                if not curv.is_real or curv == sp.sympify("+oo") or curv == sp.sympify("zoo") or curv == sp.sympify("nan"):
                    curv = curvMax
            Delta_t = (rangeValues[1]-rangeValues[0])/getNum(curv, curvMax, rangeValues[2], rangeValues[3])
            tNum += Delta_t
            if tNum <= rangeValues[1]:
                tRange = np.append(tRange, tNum)
            else:
                tRange = np.append(tRange, rangeValues[1])

    for i in reversed(range(1, len(tRange))):
        if tRange[i] == tRange[i-1]:
            tRange = np.delete(tRange, i)
    return tRange

def main():
    startTime = time.time()
    
    mirrorName = "Mirror" #placeholder name. Beware of only using valid string characters
    t = sp.Symbol('t', real = True) #mirror
    tPy = "t"
    xMirror_t = 4*sp.cos(t)*sp.cos(t)*sp.cos(t) #placeholder function
    yMirror_t = 4*sp.sin(t)*sp.sin(t)*sp.sin(t) #placeholder function
    tRangeValuesList = [(0, 2*np.pi, 4000, 1500)] #placeholder range and densities. The tuples are (start, stop, numMx, numMin[optional]), (extension_start, extension_stop, extension_numMax, extension_numMin[optional]) etc. Note that it has to be such that start < stop, etc
    tRange = generateRange(tRangeValuesList, True, xMirror_t, yMirror_t, tPy, t)
    #tRangePlot = np.linspace(0, 2*np.pi, num=100) #full parameter range to have a smooth plot of the curve, albeit doing the reflection calculations only for the limited interval of tRange
    tRangePlot = tRange

    toBeMirroredName = "ToBeMirrored" #Placeholder name. Beware of only using valid string characters
    q = sp.Symbol('q', real = True) #to be mirrored
    qPy = "q"
    xToBeMirrored_q = 4*sp.cos(q) #placeholder function
    yToBeMirrored_q = 4*sp.sin(q) #placeholder function
    qRangeValuesList = [(0, 2*np.pi, 200)] #increasing these num values vastly increases computation time
    #qIntervals = (sp.Interval(firstPieceStart, firstPieceStop), sp.Interval.Lopen(secondPieceStart, secondPieceStop), etc)
    qRange = generateRange(qRangeValuesList)

    plt.figure(num=0, dpi=150)

    xMirrorList, yMirrorList = points(xMirror_t, yMirror_t, tPy, t, tRangePlot)
    plt.plot(xMirrorList, yMirrorList, '.') if (isinstance(xMirror_t, sp.Piecewise) or isinstance(yMirror_t, sp.Piecewise)) else plt.plot(xMirrorList, yMirrorList)
    xToBeMirroredList, yToBeMirroredList = points(xToBeMirrored_q, yToBeMirrored_q, qPy, q, qRange)
    plt.plot(xToBeMirroredList, yToBeMirroredList, '.') if (isinstance(xToBeMirrored_q, sp.Piecewise) or isinstance(yToBeMirrored_q, sp.Piecewise)) else plt.plot(xToBeMirroredList, yToBeMirroredList)

    manager = multiprocessing.Manager()
    mirroredShared = manager.list()

    xSegmentList, ySegmentList, qNum_1List, qNum_2List = segment(xToBeMirrored_q, yToBeMirrored_q, qPy, q, qRange) #, qIntervals)
    nProcesses, activeProcesses = 10, []
    for i in range(0, nProcesses):
        process = multiprocessing.Process(target=mirror, args=(xSegmentList, ySegmentList, qNum_1List, qNum_2List, xMirror_t, yMirror_t, xToBeMirrored_q, yToBeMirrored_q, t, tRange, tPy, qPy, i, nProcesses, mirroredShared))
        activeProcesses.append(process)
        process.start()
    for process in activeProcesses:
        process.join()
    xMirroredList, yMirroredList = [], []
    for mirroredPoints in mirroredShared:
        xMirroredList.append(mirroredPoints[0])
        yMirroredList.append(mirroredPoints[1])
    plt.plot(xMirroredList, yMirroredList, '.')

    with open(f'{toBeMirroredName}_from_{mirrorName}.csv', "w+") as file1:
        for xMirrored, yMirrored in zip(xMirroredList, yMirroredList):
            file1.write(f"{str(xMirrored)},{str(yMirrored)}{chr(10)}")

    plt.xlim(-8, 8) #placeholder values
    plt.ylim(-8, 8) #placeholder values
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f'{toBeMirroredName}_from_{mirrorName}.png', dpi=600)
    endTime = time.time()
    executionTime = endTime - startTime
    print(f"Execution time: {executionTime} seconds.")
    plt.show()


if __name__ == '__main__':
    main()
