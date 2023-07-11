import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from itertools import islice
import multiprocessing
import time


def notReal(x_t, tPy, t, tNum):
    t = sp.Symbol('t', real = True) if tPy == "t" else sp.Symbol('q', real = True) 
    if isinstance(x_t, sp.Basic):
        if not x_t.subs(t, tNum).is_real:
            return True
    return False

def returnValue(x_t, tPy, t, tNum):
    t = sp.Symbol('t', real = True) if tPy == "t" else sp.Symbol('q', real = True)
    if isinstance(x_t, sp.Basic):
        return x_t.subs(t, tNum)
    return x_t

def returnValueNoAss(x_t, t, tNum):
    if isinstance(x_t, sp.Basic):
        return x_t.subs(t, tNum)
    return x_t

def curve(x_t, y_t, tPy, t, tNum):
    if notReal(x_t, tPy, t, tNum) or notReal(y_t, tPy, t, tNum):
        return 0, 0, False
    return returnValue(x_t, tPy, t, tNum), returnValue(y_t, tPy, t, tNum), True

def getCoeffPerp(x_t, y_t, xPrime_t1, yPrime_t1, tPy, t, t1, tNum):
    yPrime = returnValueNoAss(yPrime_t1, t1, tNum)
    xPrime = returnValueNoAss(xPrime_t1, t1, tNum)
    if not (yPrime == sp.sympify("nan") or yPrime == sp.sympify("+oo") or yPrime == sp.sympify("-oo")): #derivative at tNum is a real number. tNun is always real so no chanche of having a derivative with non-zero immaginary part
        yPrime_t = sp.diff(y_t, t) if isinstance(y_t, sp.Basic) else 0 #calculate the derivative with the real=True assumption, such that it returns just a number
        yPrime = returnValue(yPrime_t, tPy, t, tNum)
    if not (xPrime == sp.sympify("nan") or xPrime == sp.sympify("+oo") or xPrime == sp.sympify("-oo")):
        xPrime_t = sp.diff(x_t, t) if isinstance(x_t, sp.Basic) else 0
        xPrime = returnValue(xPrime_t, tPy, t, tNum)
    
    if yPrime == sp.sympify("nan") or xPrime == sp.sympify("nan"):
        return np.nan
    elif yPrime != 0 and xPrime != 0: #they can't be both inf so no worries about that
        return -xPrime/yPrime #-1/coeff tan
    elif xPrime == 0:
        return 0
    else: #yPrime == 0
        return np.inf

def segment(xToBeMirrored_t, yToBeMirrored_t, tPy, t, tRange):
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
    elif abs(coeff - abs((yToBeMirrored_2 - yToBeMirrored_1))/abs((xToBeMirrored_2 - xToBeMirrored_1))) <= 0.00001 and not linIndip(xMirror, yMirror, xToBeMirrored_1, yToBeMirrored_1, xToBeMirrored_2, yToBeMirrored_2):
        return (True, (xToBeMirrored_1, yToBeMirrored_1), (xToBeMirrored_2, yToBeMirrored_2))
    return (False, (0, 0), (0, 0))

def mirror(xSegmentList, ySegmentList, qNum_1List, qNum_2List, xMirror_t, yMirror_t, xToBeMirrored_q, yToBeMirrored_q, t, tRange, tPy, qPy, currentProcess, nProcesses, mirroredShared):
    tRange = islice(tRange, currentProcess, len(tRange), nProcesses)
    t1 = sp.Symbol('t')
    xMirror_t1, yMirror_t1 = xMirror_t.subs(t, t1), yMirror_t.subs(t, t1)
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
            mirroredShared.append(calcSimm(xMirror, yMirror, i[0], i[1]))

def calcSimm(xMirror, yMirror, xToBeMirroredIntersection, yToBeMirroredIntersection):
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

def generateRange(rangeValuesList):
    return np.hstack([np.linspace(start, stop, num=num) for start, stop, num in rangeValuesList])

def main():
    startTime = time.time()
    
    mirrorName = "Mirror" #placeholder name. Beware of only usig valid string characters
    t = sp.Symbol('t', real = True) #mirror
    tPy = "t"
    xMirror_t = 4*sp.cos(t)*sp.cos(t)*sp.cos(t) #placeholder function
    yMirror_t = 4*sp.sin(t)*sp.sin(t)*sp.sin(t) #placeholder function
    tRangeValuesList = [(0, 2*np.pi, 2500)] #placeholder range and density. The tuples are (start, stop, num), (extension_start, extension_stop, extension_num) etc
    tRange = generateRange(tRangeValuesList)
    #tRangePlot = np.linspace(0, 2*np.pi, num=100) #full parameter range to have a smooth plot of the curve, albeit doing the reflection calculations only for the limited interval of tRange
    tRangePlot = tRange

    toBeMirroredName = "ToBeMirrored" #Placeholder name. Beware of only using valid string characters
    q = sp.Symbol('q', real = True) #to be mirrored
    qPy = "q"
    xToBeMirrored_q = 4*sp.cos(q) #placeholder function
    yToBeMirrored_q = 4*sp.sin(q) #placeholder function
    qRangeValuesList = [(0, 2*np.pi, 200)] #increasing these num values vastly increases computation time
    qRange = generateRange(qRangeValuesList)

    plt.figure(num=0, dpi=150)

    xMirrorList, yMirrorList = points(xMirror_t, yMirror_t, tPy, t, tRangePlot)
    plt.plot(xMirrorList, yMirrorList)
    xToBeMirroredList, yToBeMirroredList = points(xToBeMirrored_q, yToBeMirrored_q, qPy, q, qRange)
    plt.plot(xToBeMirroredList, yToBeMirroredList)

    manager = multiprocessing.Manager()
    mirroredShared = manager.list()

    xSegmentList, ySegmentList, qNum_1List, qNum_2List = segment(xToBeMirrored_q, yToBeMirrored_q, qPy, q, qRange)
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
