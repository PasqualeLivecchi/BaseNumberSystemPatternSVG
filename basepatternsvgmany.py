from flask import Markup
import numpy, time, concurrent.futures, os, sys, traceback
from itertools import chain
from collections import OrderedDict,namedtuple
from functools import lru_cache

MINBASE = 10
MAXBASE = MINBASE+1
COEFFICIENT = 2
OP = "*"
# OPLIST = ["*","+"]


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
    ratioview2text = max(width,height)//15
    svgtext = [f"<text fill='#ae76e4' font-size='{ratioview2text}' x='{xypos[0]}' y='{xypos[1]}' transform='rotate(90 {xypos[0]},{xypos[1]})'>{numkey}</text>" for numkey, xypos in ptlocations.items()]
    htmlsvglst = [Markup(f"<div class='patterndiv'><div id='equationlabel{base}{operation}{coefficient}' class='equationlabel'>Base:{base} Operation:{operation} Coefficient: {coefficient} Equation: {base}{operation}{coefficient} </div><svg class='svg-content' width='{width}' height='{height}' viewBox='-{ratioview2text},-{ratioview2text},{width+ratioview2text*2},{height+ratioview2text*2}'>{''.join(svgtext)}<polygon points='{' '.join(ptl)}' fill='none' stroke='black' stroke-width='2'></polygon>")]
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
        abpstr = f"{abpstr}`{single_op_base_patterns(i,*bpi)}"
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
    with concurrent.futures.ThreadPoolExecutor() as threxecutor:
        basepatterninfo = namedtuple('bpi',['base','operation','coefficient','width','height','rotation','increment'])
        html = "<!doctype html><html lang='en'>"
        html += """<head>
            <meta charset="utf-8">
            <meta name="viewport" content="width=device-width, initial-scale=1">
        <style>
            .patterndiv {position: relative;}
            .equationlabel {top:40px; text-align: center;}
            .svg-content {margin-top: 8px; margin-left: 15px; padding: 6px; animation: spin 5s linear infinite; }
            .polygon-patterns {stroke-dasharray: 50; stroke-width: 2; animation: draw 5s linear reverse infinite;}
            @keyframes draw {from {stroke-dashoffset: 0;} to {stroke-dashoffset: 1000;}}
            @keyframes spin {from {transform: rotate(0deg);} to {transform: rotate(360deg);}}
        </style>
        </head>
        <body>
            <div id="topbar">
                <form id="symbolpatternform" method="POST" class="form">
                    <input id="csrf_token" name="csrf_token" type="hidden" value="ImQ3NDBiMTAwZjk1ZmFlZDg4ZjdlNWUwMGFhMjg1Y2NmODI3MTY1ZjYi.YYxx7Q.-r827NgLdpZtxaFW9Knx5XEco_s"> 
                    <label class="form-control-label" for="base">Base Number System</label>
                    <input class="form-control" id="formbaseid" name="base" type="text" value="10">
                    <label class="form-control-label" for="operation"> Operation</label>
                    <select class="form-control" id="formopid" name="operation">
                        <option selected value="*"> Multiplication *</option>
                        <option value="//">Division /</option>
                        <option value="+">Addition +</option>
                        <option value="-">Subtraction -</option>
                    </select>
                    <label class="form-control-label" for="coefficient"> Coefficient</label>
                    <input class="form-control" id="formcoefid" name="coefficient" type="text" value="3">
                    <button id="togglerotation" type="button" name="rotation" value="Rotate">Rotate</button>
                    <button id="toggleincrement" type="button" name="increment" value="Increment">Increment</button>
                </form>
            </div>"""
        htmllst,baserange = [],range(max(MINBASE,3),MAXBASE+1)
        def creatsvghtmls(base):
            try:
                bpi = basepatterninfo(base,OP,COEFFICIENT,666,666,False,False)
                htmllst.extend(create_svg_html(*bpi))
            except Exception:
                error = sys.exc_info()[0](traceback.format_exc())
                print(f"error with b op c: {base} {op} {COEFFICIENT} {error}")
                return
        try:
            threxecutor.map(creatsvghtmls,baserange)
        except Exception:
            error = sys.exc_info()[0](traceback.format_exc())
            print(f"error during threxecution {error}")
        html += "".join(htmllst)
        html += """\n<script type="text/javascript">
                    function refreshsvg(el) {
                        var svginviewlist = [];
                        var svgs = document.querySelectorAll('.polygon-patterns')
                        console.log("len: " + svgs.length.toString());
                        for (var i=0; i<svgs.length;i++){
                            let rect = svgs[i].getBoundingClientRect();
                            if (rect.top >= -911 && rect.left >= -666 && rect.bottom <= ((window.innerHeight || document.documentElement.clientHeight)+911) && rect.right <= ((window.innerWidth || document.documentElement.clientWidth)+666)){
                                svgs[i].style.display = "block";
                                svginviewlist.push(svgs[i]);
                            } else {
                                svgs[i].style.display = "none";
                                svgs[i].parentElement.style.animation = "";
                            }
                        }
                        return svginviewlist
                    }

                    function svgrotate(el){
                        var inviewlist = refreshsvg(document);
                        console.log("svgrotate"+inviewlist.length);
                        for (var i=0; i<inviewlist.length; i++){
                            let rotateaniisoff = inviewlist[i].parentElement.style.animation === "";
                            console.log("ani before" + rotateaniisoff);
                            if (rotateaniisoff){
                                inviewlist[i].parentElement.style.animation = "spin 4s linear infinite";
                            }else{
                                inviewlist[i].parentElement.style.animation = "";
                            }
                        }
                    }

                    // document.addEventListener('scroll', refreshsvg);
                    document.addEventListener('load', refreshsvg);
                    document.addEventListener('resize', refreshsvg);
                    document.getElementById('togglerotation').addEventListener('click', svgrotate);
                </script></body></html>"""
    return html

if __name__ == '__main__':
    # h = html()
    t0 = time.time()
    h = svghtml()
    print(f"svghtml:{time.time()-t0}")
    # print(h)
    with open("svgmodelmany.html",'w') as w:
        w.write(h)