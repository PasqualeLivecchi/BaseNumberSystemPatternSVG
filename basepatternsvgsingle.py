# Script that creates a html file that can be opened in a browser to show an svg image containing all patterns in the 
# base number system, operation, and coefficient equation. See for an explanation of how the patterns are found:
# https://youtu.be/kVD-6LG0-h0
from flask import Markup
import numpy, time, concurrent.futures, os, sys, traceback
from itertools import chain
from collections import OrderedDict,namedtuple
from functools import lru_cache



# BASE_NUMBER_SYSTEM number works best when using numbers that range from a min number of 2 to max number of around 600. 
# Much higher numbers can be used but over 100 or so then the numbers surrounding the graphic become unreadable,
# over 1000 and the script will take up to 3 seconds to finish. Higher than 2000 and it will take 20-30 seconds for the script to finish.
BASE_NUMBER_SYSTEM = 66
# OPERATION is the type of mathematic operation to find patterns for possible values are multiplication '*', division '/',
# addition '+', subtraction '-', AND '&', OR '|', XOR '^', fourier 'f', bitshift left '<<', bitshift right '>>', 
# some do not work and will give errors, its best to stick to multiplication '*'.
OPERATION = '*'
# COEFFICIENT is what number is being continously multiplied, 
# it can be any number from 1-9
COEFFICIENT = 3 


# USERS DON"T NEED TO CHANGE ANYTHING PAST THIS POINT
basepatterninfo = namedtuple('bpi',['base','operation','coefficient','width','height','rotation','increment'])
htmllst = []

class IndexedSet:
    def __init__(self, other=None):
        self.item_index_map = dict()
        self.item_list = []
        if other:
            self.update(other)

    def add(self, item):
        if item not in self.item_index_map:
            self.item_index_map[item] = len(self.item_list)
            self.item_list.append(item)

    def update(self, *others):
        if not others:
            return
        elif len(others) == 1:
            other = others[0]
        else:
            other = chain(others)
        for o in other:
            self.add(o)

    def __iter__(self):
        return (item for item in self.item_list if item is not None)


def create_svg_html(*bpi):
    base,operation,coefficient,width,height,hasspiral,increment = bpi
    abp = all_base_patterns(*bpi)
    ptlocations = point_locations(*bpi)
    pps = pattern_paths(abp,ptlocations)
    htmlsvglst,i,colormkuplst = [],1,[]
    while len(pps) > 0:
        od = pps.pop()
        svgpts = []
        odv = od.values()
        for pt in odv:
            if len(odv) == 1:
                percent2w,percent2h = width//50, height//50
                # this should make a square around a point when there is a repeating pattern
                svgpts.append(f"{pt[0]-percent2w},{pt[1]-percent2h}")
                svgpts.append(f"{pt[0]+percent2w},{pt[1]-percent2h}")
                svgpts.append(f"{pt[0]+percent2w},{pt[1]+percent2h}")
                svgpts.append(f"{pt[0]-percent2w},{pt[1]+percent2h}")
                svgpts.append(f"{pt[0]-percent2w},{pt[1]-percent2h}")
            else:
                svgpts.append(f"{pt[0]},{pt[1]}")
        color = str(hex(time.time_ns()*i))[-6:]
        mkupstr = f"<polygon class='polygon-patterns' points='{' '.join(svgpts)}' fill='none' stroke='#{color}' stroke-width='2'></polygon>"
        i += 1
        colormkuplst.append(Markup(f"{mkupstr}"))
    ptl = [f"{pt[0]},{pt[1]}" for pt in ptlocations.values()]
    ratioview2text = max(width,height)//35
    svgtextlst = []
    # svgtextlst displays the number in its appropriate position around the circle.
    svgtextlst = [f"<text fill='#ae76e4' font-size='{ratioview2text}' x='{xypos[0]}' y='{xypos[1]}' transform='rotate(90 {xypos[0]},{xypos[1]})'>{numkey}</text>" for numkey, xypos in ptlocations.items()]
    htmlsvglst = [Markup(f"<div class='patterndiv'><div id='equationlabel{base}{operation}{coefficient}' class='equationlabel'>Base Number System: {base} Operation: {operation} Coefficient: {coefficient} Equation: {base} {operation} {coefficient} Each color is a repeating pattern that occurs when the operation {operation} {coefficient} is repeated and then the result is reduced by addition into a single digit (number less than the base number).</div><svg class='svg-content' width='{width}' height='{height}' viewBox='-{ratioview2text},-{ratioview2text},{width+ratioview2text*2},{height+ratioview2text*2}'>{''.join(svgtextlst)}<polygon points='{' '.join(ptl)}' fill='none' stroke='black' stroke-width='2'></polygon>")]
    htmlsvglst.extend(colormkuplst)
    htmlsvglst.append(Markup("</svg></div>"))
    return htmlsvglst

def pattern_paths(abp,ptlocations):
    plst = []
    dloc = ptlocations
    while len(abp) > 0:
        pstr = abp.pop()
        l = [int(num) for num in pstr.split(',')]
        d = OrderedDict.fromkeys(l)
        haszero = False
        for k in d.keys():
            # haszero = k == 0
            # if haszero:
            #     print('haszero')
            #     break
            d[k] = dloc[k]
        if not haszero:
            plst.append(d)
    return plst

@lru_cache(maxsize=10000)
def point_locations(*bpi):
    base,operation,coefficient,w,h,hasspiral,increment = bpi
    ptdict = {}
    angle = 360/(base-1)
    hw = w//2
    hh = h//2
    dist = (w+h)//4
    spiral = max(dist//(base-1),1)
    for i in range(1,base):
        ang = angle*i#(angle*i+180) % 360 - 180
        ptdict[i] = ((dist * numpy.cos(numpy.radians(ang))+hw),(dist * numpy.sin(numpy.radians(ang))+hh))
        dist = dist - spiral if hasspiral else dist
    return ptdict

@lru_cache(maxsize=10000)
def d2b(d,base):
    quot,rem = d//base,d%base
    d2bstr = f'`{rem}'
    while quot > 0:
        rem,quot=quot%base,quot//base
        d2bstr = f'{d2bstr}`{rem}'
    return f'{base}:{d2bstr}'

@lru_cache(maxsize=10000)
def digitalroot(bdstr):
    i,b,dstr,dstrmax,bdsplit,dsum = 0,0,"",0,[],0
    bdsplit = bdstr.split(":")
    b,dstr,dstrmax = int(bdsplit[0]),bdsplit[1],dstr.count("`")
    while i <= dstrmax:
        dsum = sum(int(d) for d in dstr.split("`")[1:])
        i += 1
    while dsum >= b:
        dsum = digitalroot(d2b(dsum,b))
    return dsum

def operationswitch(digroot,op,coef):
    return {
        'mult': digroot*coef,
        'multiply': digroot*coef,
        '*': digroot*coef,
        'div': digroot//coef,
        'divide': digroot//coef,
        '/': digroot//coef,
        'add': digroot+coef,
        'plus': digroot+coef,
        '+': digroot+coef,
        'sub': digroot-coef,
        'minus': digroot-coef,
        '-': digroot-coef,
        'and': digroot&coef,
        '&': digroot&coef,
        'xor': digroot^coef,
        '^': digroot^coef,
        'or': digroot|coef,
        '|': digroot|coef,
        'fourier': int(digroot*2**(coef*4)),
        'f': int(digroot*2**(coef*4)),
        'bitshiftright': digroot>>coef,
        '>>': digroot>>coef,
        'bitshiftleft': digroot<<coef,
        '<<': digroot<<coef
    }[op]

@lru_cache(maxsize=10000)
def single_op_base_patterns(start,*bpi):
    base,op,coef,w,h,hasspiral,increment = bpi
    single_op,i,_next = "",1,start
    while i < base:
        single = digitalroot(d2b(_next,base))
        single_op = f"{single_op},{single}"
        _next = operationswitch(single,op,coef)
        i += 1
    return single_op

@lru_cache(maxsize=10000)
def all_base_patterns(*bpi):
    base,op,coef,w,h,hasspiral,increment = bpi
    i,k,abpstr,dupels = 1,0,"",[]
    while i < base:
        abpstr = f"{abpstr}`{single_op_base_patterns(i, *bpi)}"
        i += 1
    abp = abpstr.split("`")[1:]
    abpis = [IndexedSet(bp.split(',')[1:]) for bp in abp]
    allbasepats = []
    for bi in range(1,base):
        for ab in abpis:
            abps = ','.join(ab)
            if str(bi) in ab:
                if abps not in allbasepats:
                    allbasepats.append(abps)
                break
    if base%2 == 1 and op == '*' and coef == 2:
        allbasepats.append(str(base-1))
    return allbasepats

# left: 21rem; top: 6rem; 
def svghtml():
    htmllst = create_svg_html(*basepatterninfo(BASE_NUMBER_SYSTEM,OPERATION,COEFFICIENT,888,888,False,False))
    html = "<!doctype html><html lang='en'>"
    html += f"""<head>
    <style>
        .patterndiv {{position: relative;}}
        .equationlabel {{top:40px; text-align: center;}}
        #extrainfodiv {{right:80%; text-align: center;}}
        #instructiondiv {{left:20%; right:20%; text-align: center;}}
        .svg-content {{position: fixed; left:24%; transform: rotate(-90deg); }}
        .polygon-patterns {{stroke-dasharray: 200; stroke-width: 6; animation: draw 12s linear reverse infinite;}}
        @keyframes draw {{from {{stroke-dashoffset: 0;}} to {{stroke-dashoffset: 5000;}}}}
        @keyframes spin {{from {{transform: rotate(0deg);}} to {{transform: rotate(360deg);}}}}
    </style>
    </head>
    <body>"""
    html += "".join(htmllst)
    return html

if __name__ == '__main__':
    # h = html()
    t0 = time.time()
    h = svghtml()
    print(f"svghtml:{time.time()-t0}")
    # print(h)
    with open("svgmodel.html",'w') as wr:
        wr.write(h)
