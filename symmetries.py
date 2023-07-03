import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from itertools import islice
import multiprocessing


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

def curva(x_t, y_t, tPy, t, tNum):
    if notReal(x_t, tPy, t, tNum) or notReal(y_t, tPy, t, tNum):
        return 0, 0, False
    return returnValue(x_t, tPy, t, tNum), returnValue(y_t, tPy, t, tNum), True

def getCoeffPerp(xPrimo_t, yPrimo_t, tPy, t, tNum):
    yPrimo = returnValue(yPrimo_t, tPy, t, tNum)
    xPrimo = returnValue(xPrimo_t, tPy, t, tNum)
    if yPrimo != 0 and xPrimo != 0:
        return -xPrimo/yPrimo #-1/coeff tan
    elif xPrimo == 0:
        return 0
    else:
        return np.inf

def segmenta(xSpecchianda_t, ySpecchianda_t, tPy, t, tRange):
    xSegmentoList = []
    ySegmentoList = []
    tNum_1List = []
    tNum_2List = []
    q = sp.Symbol('q', real = True)
    for i in range(0, len(tRange)-1):
        tNum_1 = tRange[i]
        tNum_2 = tRange[i+1]
        xSpecchianda_1, ySpecchianda_1, reale_1 = curva(xSpecchianda_t, ySpecchianda_t, tPy, t, tNum_1)
        xSpecchianda_2, ySpecchianda_2, reale_2 = curva(xSpecchianda_t, ySpecchianda_t, tPy, t, tNum_2)
        if not reale_1 or not reale_2:
            continue
        xSegmentoList += [xSpecchianda_1+q*(xSpecchianda_2-xSpecchianda_1)]
        ySegmentoList += [ySpecchianda_1+q*(ySpecchianda_2-ySpecchianda_1)]
        tNum_1List += [tNum_1]
        tNum_2List += [tNum_2]
    return xSegmentoList, ySegmentoList, tNum_1List, tNum_2List

def intersez(xSpecchio, ySpecchio, coeff, xSegmentoList, ySegmentoList, qNum_1List, qNum_2List, tNum, xSpecchio_t, ySpecchio_t, xSpecchianda_q, ySpecchianda_q, tPy, qPy):
    intersezioni = []
    q = sp.Symbol('q', real = True)
    t = sp.Symbol('t', real = True)
    r = sp.Symbol('r', real = True)
    for i in range(0, len(xSegmentoList)-1):
        xRetta_r = xSpecchio+r if coeff != np.inf else xSpecchio
        yRetta_r = ySpecchio+r*coeff if coeff != np.inf else r

        intersezione = None
        try: #could be done more efficiently with some if/else
            intersezione = sp.linsolve([xSegmentoList[i]-xRetta_r, ySegmentoList[i]-yRetta_r], q, r)
            qNum, rNum = list(intersezione)[0] #if the same, qNum is symbolical. If without intersections, length is 0
            if r in qNum.free_symbols or q in qNum.free_symbols:
                coin, specchianda_1Tupla, specchianda_2Tupla = coincidenti(xSpecchio_t, ySpecchio_t, t, tNum, xSpecchianda_q, ySpecchianda_q, q, qNum_1List[i], qNum_2List[i], coeff, tPy, qPy)
                if coin:
                    intersezioni += [specchianda_1Tupla]
                    intersezioni += [specchianda_2Tupla]
            elif 0 <= qNum <= 1:
                intersezioni += [(xSegmentoList[i].subs(q, qNum), ySegmentoList[i].subs(q, qNum))]
        except ValueError:
            continue
        except IndexError:
            #probably just a continue would do
            if len(list(intersezione)) != 0:
                coin, specchianda_1Tupla, specchianda_2Tupla = coincidenti(xSpecchio_t, ySpecchio_t, t, tNum, xSpecchianda_q, ySpecchianda_q, q, qNum_1List[i], qNum_2List[i], coeff, tPy, qPy)
                if coin:
                    intersezioni += [specchianda_1Tupla]
                    intersezioni += [specchianda_2Tupla]
    return intersezioni

def linIndip(xSpecchio, ySpecchio, xSpecchianda_1, ySpecchianda_1, xSpecchianda_2, ySpecchianda_2):
    a = sp.Symbol('a', real = True)
    b = sp.Symbol('b', real = True)
    dipendenza = list(sp.linsolve([a*(xSpecchianda_1-xSpecchio)+b*(xSpecchianda_2-xSpecchio), a*(ySpecchianda_1-ySpecchio)+b*(ySpecchianda_2-ySpecchio)], a, b))
    aNum, bNum = dipendenza[0]
    if isinstance(aNum, int) and isinstance(bNum, int): #maybe useless
        if aNum == 0 and bNum == 0:
            return True
    return False

def coincidenti(xSpecchio_t, ySpecchio_t, t, tNum, xSpecchianda_q, ySpecchianda_q, q, qNum_1, qNum_2, coeff, tPy, qPy):
    xSpecchio = returnValue(xSpecchio_t, tPy, t, tNum)
    ySpecchio = returnValue(ySpecchio_t, tPy, t, tNum)
    xSpecchianda_1 = returnValue(xSpecchianda_q, qPy, q, qNum_1)
    xSpecchianda_2 = returnValue(xSpecchianda_q, qPy, q, qNum_2)
    ySpecchianda_1 = returnValue(ySpecchianda_q, qPy, q, qNum_1)
    ySpecchianda_2 = returnValue(ySpecchianda_q, qPy, q, qNum_2)
    if xSpecchianda_1 == xSpecchianda_2:
        if not linIndip(xSpecchio, ySpecchio, xSpecchianda_1, ySpecchianda_1, xSpecchianda_2, ySpecchianda_2):
            return (True, (xSpecchianda_1, ySpecchianda_1), (xSpecchianda_2, ySpecchianda_2))
    elif abs(coeff - abs((ySpecchianda_2 - ySpecchianda_1))/abs((xSpecchianda_2 - xSpecchianda_1))) <= 0.0001 and not linIndip(xSpecchio, ySpecchio, xSpecchianda_1, ySpecchianda_1, xSpecchianda_2, ySpecchianda_2):
        return (True, (xSpecchianda_1, ySpecchianda_1), (xSpecchianda_2, ySpecchianda_2))
    return (False, (0, 0), (0, 0))

def specchia(xSegmentoList, ySegmentoList, qNum_1List, qNum_2List, xSpecchio_t, ySpecchio_t, xSpecchianda_q, ySpecchianda_q, t, tRange, tPy, qPy, currentProcess, nProcesses, specchiataShared):
    tRange = islice(tRange, currentProcess, len(tRange), nProcesses)
    xPrimo_t = sp.diff(xSpecchio_t, t) if isinstance(xSpecchio_t, sp.Basic) else 0
    yPrimo_t = sp.diff(ySpecchio_t, t) if isinstance(ySpecchio_t, sp.Basic) else 0
    for tNum in tRange:
        print(tNum)
        xSpecchio, ySpecchio, reale = curva(xSpecchio_t, ySpecchio_t, tPy, t, tNum)
        if not reale:
            continue
        coeff = getCoeffPerp(xPrimo_t, yPrimo_t, tPy, t, tNum)
        intersezioni = intersez(xSpecchio, ySpecchio, coeff, xSegmentoList, ySegmentoList, qNum_1List, qNum_2List, tNum, xSpecchio_t, ySpecchio_t, xSpecchianda_q, ySpecchianda_q, tPy, qPy)
        for i in intersezioni:
            specchiataShared.append(calcolaSimm(xSpecchio, ySpecchio, i[0], i[1]))

def calcolaSimm(xSpecchio, ySpecchio, xSpecchiandaIntersezione, ySpecchiandaIntersezione):
    return (2*xSpecchio-xSpecchiandaIntersezione, 2*ySpecchio-ySpecchiandaIntersezione)

def punti(x_t, y_t, tPy, t, tRange):
    x_tList = []
    y_tList = []
    for tNum in tRange:
        if notReal(x_t, tPy, t, tNum) or notReal(y_t, tPy, t, tNum):
            continue
        x_tList += [returnValue(x_t, tPy, t, tNum)]
        y_tList += [returnValue(y_t, tPy, t, tNum)]
    return x_tList, y_tList

def main():
    specchioNome = "Specchio"
    t = sp.Symbol('t', real = True) #mirror
    tPy = "t"
    xSpecchio_t = 4*sp.cos(t)*sp.cos(t)*sp.cos(t)
    ySpecchio_t = 4*sp.sin(t)*sp.sin(t)*sp.sin(t) #placeholder function
    tRange = np.linspace(0, 2*np.pi, num=2500) #placeholder function
    #tRange.extend(np.linspace())

    specchiandaNome = "Specchianda"
    q = sp.Symbol('q', real = True) #to be mirrored
    qPy = "q"
    xSpecchianda_q = 4*sp.cos(q) #placeholder function
    ySpecchianda_q = 4*sp.sin(q) #placeholder function
    qRange = np.linspace(0, 2*np.pi, num=100)
    #qRange.extend(np.linspace())

    plt.figure(num=0, dpi=150)

    xSpecchioList, ySpecchioList = punti(xSpecchio_t, ySpecchio_t, tPy, t, tRange)
    plt.plot(xSpecchioList, ySpecchioList)
    xSpecchiandaList, ySpecchiandaList = punti(xSpecchianda_q, ySpecchianda_q, qPy, q, qRange)
    plt.plot(xSpecchiandaList, ySpecchiandaList)

    manager = multiprocessing.Manager()
    specchiataShared = manager.list()

    xSegmentoList, ySegmentoList, qNum_1List, qNum_2List = segmenta(xSpecchianda_q, ySpecchianda_q, qPy, q, qRange)
    nProcesses, activeProcesses = 12, [] #for optimum performance don't exceed the number of physical cores of your CPU (x2 if SMT (AMD) or Hyper-Threading (Intel) is enabled)
    for i in range(0, nProcesses):
        process = multiprocessing.Process(target=specchia, args=(xSegmentoList, ySegmentoList, qNum_1List, qNum_2List, xSpecchio_t, ySpecchio_t, xSpecchianda_q, ySpecchianda_q, t, tRange, tPy, qPy, i, nProcesses, specchiataShared))
        activeProcesses.append(process)
        process.start()
    for process in activeProcesses:
        process.join()
    xSpecchiataList, ySpecchiataList = [], []
    for specchiataPoints in specchiataShared:
        xSpecchiataList.append(specchiataPoints[0])
        ySpecchiataList.append(specchiataPoints[1])
    plt.plot(xSpecchiataList, ySpecchiataList, '.')

    with open(f'{specchiandaNome}_from_{specchioNome}.csv', "w+") as file1:
        for xSpecchiata, ySpecchiata in zip(xSpecchiataList, ySpecchiataList):
            file1.write(f"{str(xSpecchiata)},{str(ySpecchiata)}{chr(10)}")

    plt.xlim(-8, 8)
    plt.ylim(-8, 8)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(f'{specchiandaNome}_from_{specchioNome}.png')
    plt.show()


if __name__ == '__main__':
    main()
