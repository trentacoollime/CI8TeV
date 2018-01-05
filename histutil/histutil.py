#-----------------------------------------------------------------------------
# Some ROOT histogram utilities
# Created: sometime in the early 21st century, HBP
# Updated: 28-May-2016 HBP - fix Table (for RGS paper)
#          09-Apr-2017 HBP - add ranking to BDT
#          29-Aug-2017 HBP - change to import ROOT
#-----------------------------------------------------------------------------
import os, sys, re
import ROOT
from glob import glob
from array import array
from string import split, strip, atoi, atof, replace, joinfields
from math import *
#-----------------------------------------------------------------------------
# Hack to suppress harmless warning
# see: https://root.cern.ch/phpBB3/viewtopic.php?f=14&t=17682
import warnings
warnings.filterwarnings( action='ignore', category=RuntimeWarning,
                         message='.*class stack<RooAbsArg\*,deque<RooAbsArg\*> >' )
warnings.filterwarnings( action='ignore', category=RuntimeWarning,
                         message='Error in <TClass::BuildRealData>: '\
                             'Cannot find any ShowMembers function for G__CINT_IOFLAGS!')
#---------------------------------------------------------------------------
def nameonly(s):
    import posixpath
    return posixpath.splitext(posixpath.split(s)[1])[0]
#---------------------------------------------------------------------------
RED    ="\x1b[0;31;48m"
GREEN  ="\x1b[0;32;48m"
YELLOW ="\x1b[0;33;48m"
BLUE   ="\x1b[0;34;48m"
MAGENTA="\x1b[0;35;48m"
CYAN   ="\x1b[0;36;48m"

BOLDRED    ="\x1b[1;31;48m"
BOLDGREEN  ="\x1b[1;32;48m"
BOLDYELLOW ="\x1b[1;33;48m"
BOLDBLUE   ="\x1b[1;34;48m"
BOLDMAGENTA="\x1b[1;35;48m"
BOLDCYAN   ="\x1b[1;36;48m"

RESETCOLOR ="\x1b[0m"    # reset to default foreground color             
#------------------------------------------------------------------------------
WIDTH  = 500    # Default canvas width
HEIGHT = 500    # Default canvas height
NDIVX  = 510

TEXTFONT = 42
# Font code = 10 * fontnumber + precision with fontnumber 1..14
# fontnumber= 1..14
# precision = 0 fast hardware fonts
# precision = 1 scalable and rotatable hardware fonts

# precision = 2 scalable and rotatable hardware fonts
# precision = 3 scalable and rotatable hardware fonts. Text size in pixels. 

Arial             = 43
FONTSIZE          = 18 # font size in pixels

TimesNewRomanBold = 33
LABEL_FONTSIZE    = 18
FONTSCALE         = 0.0027

TimesNewRoman2    = 12
TITLE_FONTSIZE    = 0.04

# Text alignment = 10 * horizontal + vertical
# For horizontal alignment:
#   1=left adjusted, 2=centered, 3=right adjusted
# For vertical alignment:
#   1=bottom adjusted, 2=centered, 3=top adjusted
ALIGN_LEFT_BOTTOM = 11
#------------------------------------------------------------------------------
def deltaPhi(phi1, phi2):
    deltaphi = phi2 - phi1;
    if fabs(deltaphi) > pi:
        return 2 * pi - fabs(deltaphi)
    else:
        return deltaphi

def setStyle():
    style = ROOT.TStyle("histutil", "histutil")
    style.SetPalette(1)
    
    # For the canvases
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(ROOT.kWhite)
    style.SetCanvasDefH(500) #Height of canvas
    style.SetCanvasDefW(500) #Width of canvas
    style.SetCanvasDefX(0)   #Position on screen
    style.SetCanvasDefY(0)

    # For the pads
    style.SetPadBorderMode(0)
    style.SetPadColor(ROOT.kWhite)
    style.SetPadGridX(ROOT.kFALSE)
    style.SetPadGridY(ROOT.kFALSE)
    style.SetGridColor(ROOT.kGreen)
    style.SetGridStyle(3)
    style.SetGridWidth(1)

    # For the frames
    style.SetFrameBorderMode(0)
    style.SetFrameBorderSize(1)
    style.SetFrameFillColor(0)
    style.SetFrameFillStyle(0)
    style.SetFrameLineColor(1)
    style.SetFrameLineStyle(1)
    style.SetFrameLineWidth(1)

    # For the histograms
    style.SetHistLineColor(ROOT.kBlack)
    style.SetHistLineStyle(0)
    style.SetHistLineWidth(2)

    style.SetEndErrorSize(2)
    #style.SetErrorX(0.)

    style.SetMarkerSize(0.4)
    style.SetMarkerStyle(20)

    # For the fit/function:
    style.SetOptFit(1)
    style.SetFitFormat("5.4g")
    style.SetFuncColor(2)
    style.SetFuncStyle(1)
    style.SetFuncWidth(1)

    # For the date:
    style.SetOptDate(0)

    # For the statistics box:
    style.SetOptFile(0)
    style.SetOptStat("")
    # To display the mean and RMS:
    # style.SetOptStat("mr") 
    style.SetStatColor(ROOT.kWhite)
    style.SetStatFont(TEXTFONT)
    style.SetStatFontSize(0.03)
    style.SetStatTextColor(1)
    style.SetStatFormat("6.4g")
    style.SetStatBorderSize(1)
    style.SetStatH(0.2)
    style.SetStatW(0.3)

    # Margins:
    style.SetPadTopMargin(0.05)
    style.SetPadBottomMargin(0.16)
    style.SetPadLeftMargin(0.20)
    style.SetPadRightMargin(0.10)

    # For the Global title:
    style.SetOptTitle(0) 
    style.SetTitleFont(TEXTFONT)
    style.SetTitleColor(1)
    style.SetTitleTextColor(1)
    style.SetTitleFillColor(10)
    style.SetTitleFontSize(0.05)

    # For the axis titles:
    style.SetTitleColor(1, "XYZ")
    style.SetTitleFont(TEXTFONT, "XYZ")
    style.SetTitleSize(0.05, "XYZ")
    style.SetTitleXOffset(1.25)
    style.SetTitleYOffset(1.40)

    # For the axis labels:
    style.SetLabelColor(1, "XYZ")
    style.SetLabelFont(TEXTFONT, "XYZ")
    style.SetLabelOffset(0.020, "XYZ")
    style.SetLabelSize(0.05, "XYZ")

    # For the axis:
    style.SetAxisColor(1, "XYZ")
    style.SetStripDecimals(ROOT.kTRUE)
    style.SetTickLength(0.03, "XYZ")
    style.SetNdivisions(NDIVX, "XYZ")
    # To get tick marks on the opposite side of the frame
    style.SetPadTickX(1)  
    style.SetPadTickY(1)

    # Change for log plots:
    style.SetOptLogx(0)
    style.SetOptLogy(0)
    style.SetOptLogz(0)

    # Postscript options:
    style.SetPaperSize(20.,20.)
    style.cd()
#------------------------------------------------------------------------------
def expo(x, fmt="%4.2f", code="#"):
    from string import replace, split, atof, atoi
    s = "%10.3e" % x
    s = replace(s, "e", " ")
    a = split(s)
    n = atof(a[0])
    m = atoi(a[1])
    num = fmt % n
    if m != 0:
        s = "%s%stimes10^{%d}" % (num, code, m)
    else:
        s = num
    return s
#------------------------------------------------------------------
class TimeLeft:
    def __init__(self, ntoys):
        self.ntoys = ntoys
        self.swatch = ROOT.TStopwatch()
    def __del__(self):
        pass
    def __call__(self, ii):
        swatch = self.swatch
        ntoys = self.ntoys
        # elapsed time since start
        s = swatch.RealTime(); swatch.Continue()
        # time/loop
        loop = ii+1
        t = s / loop
        # time left
        s = t * (ntoys-loop)
        h = int(s / 3600) 
        s = s - 3600*h
        m = int(s / 60)
        s =  s - 60*m
        rec = "%2.2d:%2.2d:%2.2d hours" % (h, m, s)                             
        return rec
#------------------------------------------------------------------------------
class Scribe(ROOT.TLatex):
    def __init__(self, xxpos, yypos, size=TITLE_FONTSIZE, font=TEXTFONT):
        ROOT.TLatex.__init__(self)
        self.SetNDC()
        self.SetTextSize(size)
        self.SetTextFont(font)
        self.SetTextAlign(12)
        
        self.xpos = xxpos
        self.ypos = yypos
        self.linewidth = 1.5*size

    def __del__(self):
        pass

    def write(self, text, xoffset=0):
        y = self.ypos
        if y < 0: return
        self.DrawLatex(self.xpos+xoffset, y, text)
        self.ypos -= self.linewidth

    def vspace(self, f=0.5):
        y = self.ypos
        if y < 0: return
        self.DrawLatex(self.xpos, y, " ")
        self.ypos -= self.linewidth * f;
#------------------------------------------------------------------------------
def addTitle(title="CMS Preliminary       L = 30fb^{-1}    #sqrt{s}=13TeV",
             size=TITLE_FONTSIZE):
    xtext = 0.200
    ytext = 0.975
    s = Scribe(xtext, ytext, size)
    s.write(title)
    return s
#------------------------------------------------------------------------------
PERCENT = [0.0230, 0.1579, 0.5000, 0.8415, 0.9770]
def percentiles(points, percent):
    pts = map(lambda x: x, points)
    pts.sort()
    np = len(pts)
    ps = []
    for p in percent:
        x = p * np
        k = int(x)
        f = x - k
        c1= pts[k]
        if k < np-1:
            c2 = pts[k+1]
            c  = c1*(1-f) + c2*f
        else:
            c = c1
        ps.append(c)
    return ps

class PercentileCurve:

    def __init__(self, size):
        self.size = size
        self.points = map(lambda x: [], size*[0])

    def __del__(self):
        pass

    def add(self, curve):
        # check if this is a histogram
        try:
            nbins = curve.GetNbinsX()
            if nbins != self.size:
                print "*** PercentileCurve - "\
                  "ERROR*** wrong number of bins (points on curve)"
                print "nbins: %d, size: %d" % (nbins, self.size)
                return False

            for ii in xrange(nbins):
                c = curve.GetBinContent(ii+1)
                self.points[ii].append(c)
        except:
            if len(curve) != self.size:
                print "*** PercentileCurve - ERROR*** wrong number of points on curve"
                print "len(curve): %d, size: %d" % (len(curve), self.size)
                return False

            for ii, y in enumerate(curve):
                self.points[ii].append(y)
        return True

    def __call__(self, percentile):
        z = []
        for ii in xrange(self.size):
            points = self.points[ii]
            points.sort()
            np = len(points)

            x = percentile * np
            k = int(x)
            f = x - k
            c1= points[k]
            if k < np-1:
                c2 = points[k+1]
                c  = c1*(1-f) + c2*f
            else:
                c = c1				
            z.append(c)
        return z

    def curves(self, percentiles):
        from array import array
        lines = []
        for p in percentiles:
            c = array('d')
            c.fromlist(self.__call__(p))
            lines.append(c)
        return lines

class StandardCurve:

    def __init__(self, size):
        self.size = size
        self.points = map(lambda x: [], size*[0])

    def __del__(self):
        pass

    def add(self, curve):

        # check if this is a histogram
        try:
            nbins = curve.GetNbinsX()
            if nbins != self.size:
                print "*** StandardCurve - ERROR*** "\
                  "wrong number of bins (points on curve)"
                print "nbins: %d, size: %d" % (nbins, self.size)
                return False

            for ii in xrange(nbins):
                c = curve.GetBinContent(ii+1)
                self.points[ii].append(c)
        except:
            if len(curve) != self.size:
                print "*** StandardCurve - ERROR*** wrong number of points on curve"
                print "len(curve): %d, size: %d" % (len(curve), self.size)
                return False

            for ii, y in enumerate(curve):
                self.points[ii].append(y)
        return True

    def __call__(self, nsigma):
        z = []
        for ii in xrange(self.size):
            points = self.points[ii]
            c = sum(points)/len(points)
            ec= sqrt(sum(map(lambda x: (x-c)**2, points))/len(points))
            if nsigma != 0: c = c + nsigma * ec
            z.append(c)
        return z

    def curves(self, sigmas):
        from array import array
        lines = []
        for p in sigmas:
            c = array('d')
            c.fromlist(self.__call__(p))
            lines.append(c)
        return lines
#------------------------------------------------------------------------------
def getarg(args, key, d):
    if args.has_key(key):
        return args[key]
    else:
        return d
#------------------------------------------------------------------------------
def mkpline(xx, y1, y2, boundary, **args):    
    color  = getarg(args, 'color',   ROOT.kYellow)
    fstyle = getarg(args, 'fstyle',  3001)
    lwidth = getarg(args, 'lwidth',  2)

    nbins = len(xx)

    # lower curve
    x = array('d')
    if type(xx) == type([]):		
        x.fromlist(xx)
    else:
        for z in xx: x.append(z)

    y = array('d')
    if type(y1) == type([]):
        y.fromlist(y1)
    else:
        for z in y1: y.append(z)

    # upper curve

    for i in xrange(nbins):
        j = nbins - i - 1
        x.append(xx[j])
        y.append(y2[j])

    x.append(xx[0])
    y.append(y1[0])

    # clip polygon
    np  = 2*nbins
    npp = 2*np
    xc = array('d')
    xc.fromlist(npp*[0])
    yc = array('d')
    yc.fromlist(npp*[0])

    if type(boundary) != type([]):
        xmin = boundary.GetBinLowEdge(1)
        xmax = boundary.GetBinLowEdge(nbins)+boundary.GetBinWidth(nbins)
        ymin = boundary.GetMinimum()
        ymax = boundary.GetMaximum()
    else:
        xmin, xmax, ymin, ymax = boundary

    np = ROOT.gPad.ClipPolygon(np, x, y, npp, xc, yc, xmin, ymin, xmax, ymax)
    pl = ROOT.TPolyLine(np, xc, yc)
    pl.SetLineColor(color)
    pl.SetLineWidth(lwidth)
    pl.SetFillColor(color)
    pl.SetFillStyle(fstyle)

    return pl
#------------------------------------------------------------------------------
def mkhist1(hname, xtitle, ytitle, nbins, xmin, xmax, **args):
    ymin   = getarg(args, 'ymin', None)
    ymax   = getarg(args, 'ymax', None)
    color  = getarg(args, 'color',   ROOT.kBlack)
    lstyle = getarg(args, 'lstyle',  1)
    lwidth = getarg(args, 'lwidth',  1)
    ndivx  = getarg(args, 'ndivx',   505)
    ndivy  = getarg(args, 'ndivy',   510)

    h = ROOT.TH1F(hname, "", nbins, xmin, xmax)		
    h.SetLineColor(color)
    h.SetLineStyle(lstyle)
    h.SetLineWidth(lwidth)

    h.SetMarkerSize(0.8)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(20)

    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleOffset(1.2);
    h.GetXaxis().SetLimits(xmin, xmax)
    h.SetNdivisions(ndivx, "X")

    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.6)
    if ymin != None: h.SetMinimum(ymin)
    if ymax != None: h.SetMaximum(ymax)
    h.SetNdivisions(ndivy, "Y")
    return h
#------------------------------------------------------------------------------
def mkhist2(hname, xtitle, ytitle,
        nbinx, xmin, xmax,	
        nbiny, ymin, ymax, **args):
    color  = getarg(args, 'color',   ROOT.kBlack)
    mstyle = getarg(args, 'mstyle',  20)
    msize  = getarg(args, 'msize',   0.5)
    ndivx  = getarg(args, 'ndivx',   505)
    ndivy  = getarg(args, 'ndivy',   505)

    h = ROOT.TH2F(hname, "", nbinx, xmin, xmax, nbiny, ymin, ymax)
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerSize(msize)
    h.SetMarkerStyle(mstyle)

        #h.GetXaxis().CenterTitle()
    h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleOffset(1.3)
    h.SetNdivisions(ndivx, "X")

        #h.GetYaxis().CenterTitle()
    h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.3)
    h.SetNdivisions(ndivy, "Y")
    return h
def fixhist2(h, xtitle=None, ytitle=None, **args):
    color  = getarg(args, 'color',   ROOT.kBlack)
    mstyle = getarg(args, 'mstyle',  20)
    msize  = getarg(args, 'msize',   0.5)
    ndivx  = getarg(args, 'ndivx',   505)
    ndivy  = getarg(args, 'ndivy',   505)

    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerSize(msize)
    h.SetMarkerStyle(mstyle)

    if xtitle != None: h.GetXaxis().SetTitle(xtitle)
    h.GetXaxis().SetTitleOffset(1.3)
    h.SetNdivisions(ndivx, "X")

    if ytitle != None: h.GetYaxis().SetTitle(ytitle)
    h.GetYaxis().SetTitleOffset(1.6)
    h.SetNdivisions(ndivy, "Y")
#------------------------------------------------------------------------------
def mkgraph(x, y, xtitle, ytitle, xmin, xmax, **args):
    ymin   = getarg(args, 'ymin', None)
    ymax   = getarg(args, 'ymax', None)
    color  = getarg(args, 'color',   ROOT.kBlack)
    lstyle = getarg(args, 'lstyle',  1)
    lwidth = getarg(args, 'lwidth',  1)
    msize  = getarg(args, 'msize',   0.5)
    mstyle = getarg(args, 'mstyle', 20)
    ndivx  = getarg(args, 'ndivx', 505)
    ndivy  = getarg(args, 'ndivy', 505)
    name   = getarg(args, 'name', None)

    if y == None:
        g = ROOT.TGraph()
    else:
        n = len(y)
        xx = x
        if type(x) == type([]):
            xx = array('d')
            for i in range(n): xx.append(x[i])
        yy = y	
        if type(y) == type([]):
            yy = array('d')
            for i in range(n): yy.append(y[i])

        g = ROOT.TGraph(n, xx, yy)

    if name != None: g.SetName(name)

    g.SetLineColor(color)
    g.SetLineStyle(lstyle)
    g.SetLineWidth(lwidth)

    g.SetMarkerColor(color)
    g.SetMarkerSize(msize)
    g.SetMarkerStyle(mstyle);

    g.GetXaxis().SetTitle(xtitle)
    g.GetXaxis().SetTitleOffset(1.2);
    g.GetXaxis().SetLimits(xmin, xmax)
    g.GetHistogram().SetNdivisions(ndivx, "X")

    g.GetYaxis().SetTitle(ytitle)
    g.GetYaxis().SetTitleOffset(1.6)
    if ymin != None and ymax != None:
        g.GetHistogram().SetAxisRange(ymin, ymax, "Y")
    g.GetHistogram().SetNdivisions(ndivy, "Y")
    return g

def mkgraphErrors(x, y, ex, ey, xtitle, ytitle, xmin, xmax, **args):
    ymin   = getarg(args, 'ymin', None)
    ymax   = getarg(args, 'ymax', None)
    color  = getarg(args, 'color',   ROOT.kBlack)
    lstyle = getarg(args, 'lstyle',  1)
    lwidth = getarg(args, 'lwidth',  1)
    ndivx  = getarg(args, 'ndivx',   505)
    ndivy  = getarg(args, 'ndivy',   510)

    n = len(y)
    xx = x
    if type(x) == type([]):
        xx = array('d')
        for i in range(n): xx.append(x[i])
    yy = y	
    if type(y) == type([]):
        yy = array('d')
        for i in range(n): yy.append(y[i])

    exx = ex
    if type(ex) == type([]):
        exx = array('d')
        for i in range(n): exx.append(ex[i])
    eyy = ey	
    if type(ey) == type([]):
        eyy = array('d')
        for i in range(n): eyy.append(ey[i])

    g = ROOT.TGraphErrors(n, xx, yy, exx, eyy)

    g.SetLineColor(color)
    g.SetLineStyle(lstyle)
    g.SetLineWidth(lwidth)

    g.SetMarkerSize(0.8)
    g.SetMarkerColor(color)
    g.SetMarkerStyle(20);

    g.GetXaxis().SetTitle(xtitle)
    g.GetXaxis().SetTitleOffset(1.2);
    g.GetXaxis().SetLimits(xmin, xmax)
    g.GetHistogram().SetNdivisions(ndivx, "X")

    g.GetYaxis().SetTitle(ytitle)
    g.GetYaxis().SetTitleOffset(1.6)
    if ymin != None and ymax != None:
        g.GetHistogram().SetAxisRange(ymin, ymax, "Y")
    g.GetHistogram().SetNdivisions(ndivy, "Y")
    return g
#------------------------------------------------------------------------------
def mkcdf(hist, minbin=1):
    c = [0.0]*(hist.GetNbinsX()-minbin+2)
    j=0
    for ibin in xrange(minbin, hist.GetNbinsX()+1):
        c[j] = c[j-1] + hist.GetBinContent(ibin)
        j += 1
    c[j] = hist.Integral()
    return c

def mkroc(name, hsig, hbkg, lcolor=ROOT.kBlue, lwidth=2, ndivx=505, ndivy=505):
    from array import array
    csig = mkcdf(hsig)
    cbkg = mkcdf(hbkg)
    npts = len(csig)
    esig = array('d')
    ebkg = array('d')
    for i in xrange(npts):
        esig.append(1 - csig[npts-1-i])
        ebkg.append(1 - cbkg[npts-1-i])
    g = ROOT.TGraph(npts, ebkg, esig)
    g.SetName(name)
    g.SetLineColor(lcolor)
    g.SetLineWidth(lwidth)

    g.GetXaxis().SetTitle("#font[12]{#epsilon_{b}}")
    g.GetXaxis().SetLimits(0,1)

    g.GetYaxis().SetTitle("#font[12]{#epsilon_{s}}")
    g.GetHistogram().SetAxisRange(0,1, "Y");

    g.GetHistogram().SetNdivisions(ndivx, "X")
    g.GetHistogram().SetNdivisions(ndivy, "Y")
    return g

def mklegend(xx, yy, xw, yw):
    lg = ROOT.TLegend(xx, yy, xx+xw, yy+yw)
    lg.SetFillColor(ROOT.kWhite)
    lg.SetTextFont(TEXTFONT)
    lg.SetBorderSize(0)
    lg.SetShadowColor(ROOT.kWhite)
    return lg
#------------------------------------------------------------------------------
class Row:
    def __init__(self, rownumber, varmap, data):
        self.row = rownumber
        self.varmap= varmap
        self.data  = data
        self.items = map(lambda x: (x[1],x[0]), self.varmap.items())
        self.items.sort()

        # Initialize row counter
        self.col = 0
        self.maxcol = len(self.items)-1

    def __del__(self):
        pass

    def __call__(self, variable):
        if not self.varmap.has_key(variable): return None
        index = self.varmap[variable]
        if len(index) == 1:
            return self.data[index[0]]
        else:
            return self.data[index[0]:index[1]+1]
        
    def __str__(self):
        strrep = "row: %d\n" % self.row
        lname = ''
        for index, name in self.items:
            for ii in index:
                v = self.data[ii]
                strvalue = ''                
                if   type(v) == type(1.0):
                    strvalue += "%12.3f" % v
                elif type(v) == type(1):
                    strvalue += "%12d" % v                        
                elif type(v) == type(""):
                    strvalue += "%12s" % v
                strrep += "%4d %-16s %s\n" % (ii, name, strvalue)
        return strip(strrep)

    # Implement Python iterator protocol	
    def __iter__(self):
        return self

    def next(self):
        if self.col > self.maxcol:
            self.col = 0
            raise StopIteration
        else:
            index, name = self.items[self.col]
            if len(index) == 1: 
                value = self.data[index[0]]
            else:
                value = self.data[index[0]:index[1]+1]
            self.col += 1
            return (name, value)

    def __len__(self):
        return len(self.varmap)

    def __getitem__(self, key):
        if type(key) != type(1): return None
        if key < -len(self.varmap): return None
        if key > len(self.varmap)-1: return None
        return self.data[key]

def tonumber(x):
    try:
        y = atof(x)
    except:
        y = x
    return y

class Table:

    def __init__(self, filename, nrows=-1):
        try:
            myfile = open(filename, 'r')
        except:
            print "*** can't read file %s" % filename
            sys.exit(0)

        # Read header and check for array variables.
        # An array variable is identified by the syntax "name number"
        # Expand header by replicating the name appending to it an number
        
        records = myfile.readlines()
        t = split(records[0])
        
        self.varname= []
        self.header = []
        self.varmap = {}
        i = 0
        while i < len(t):
            name = t[i] # should be a name
            try:
                size = atoi(name)
                sys.exit('** wrong header syntax\n'\
                         '** found an integer in column %d where a string was '\
                         'expected\n' % i)
            except:
                pass
            
            # check if this is an array variable by looking ahead
            size = 1
            array_type = False
            if i < len(t)-1:
                try:
                    size = atoi(t[i+1])
                    array_type = True
                except:
                    pass

            if array_type:
                # this is an array variable so expand header names
                # (if size > 1),
                # but first cache starting position of array
                # and its size
                self.varmap[name] = []
                if size > 1:
                    self.varmap[name] = [len(self.header),
                                         len(self.header)+size-1]
                    for j in xrange(size):
                        newname = '%s[%d]' % (name, j)
                        self.header.append(newname)
                else:
                    self.varmap[name] = [len(self.header)]
                    self.header.append(name)
                i += 2
            else:
                # this is scalar variable so just add it to the
                # header as-is
                self.varmap[name] = [len(self.header)]
                self.header.append(name)
                i += 1
            # cache original names and size
            self.varname.append((name, size))
            
        ## # print some stuff
        ## print 
        ## lname = ''
        ## for ii, name in enumerate(self.header):
        ##     basename = split(name, '[')[0]
        ##     if basename != lname:
        ##         extra = self.varmap[basename]
        ##     else:
        ##         extra = ''
        ##     lname = basename
        ##     print '%5d\t%-16s\t%s' % (ii, self.header[ii], extra) 

        # load data into memory
        rownumber = 0
        self.data = []
        index = 1
        for record in records[1:]:
            # Convert to numbers
            record = map(tonumber, split(record))
            self.data.append(record)
            rownumber += 1
            if nrows > 0:
                if rownumber >= nrows:
                    break
        myfile.close()

        # Initialize row counter
        self.rownumber = 0
        self.maxrow = len(self.data)-1
        self.maxcol = len(self.header)-1
        
        # Create a name to index map for rows
        self.rowmap = {} # empty map
        for index in xrange(len(self.data)):
            self.rowmap['%d' % index] = index            

    def __del__(self):
        pass

    def __call__(self, rownumber, variable=None):
        if rownumber < 0: return None
        if rownumber > self.maxrow: return None
        if variable == None:
            return Row(rownumber, self.varmap, self.data[rownumber])
        else:
            if not self.varmap.has_key(variable): return None
            index = self.varmap[variable]
            if len(index) == 1:
                return self.data[rownumber][index[0]]
            else:
                return self.data[rownumber][index[0]:index[1]+1]

    # Implement Python iterator protocol
    def __iter__(self):
        return self

    def next(self):
        if self.rownumber > self.maxrow:
            self.rownumber = 0
            raise StopIteration
        else:
            data = Row(self.rownumber,
                       self.varmap,
                       self.data[self.rownumber])
            self.rownumber += 1
            return data

    def row(self, rownumber):
        if rownumber > self.maxrow: return None
        data = Row(rownumber, self.varmap, self.data[rownumber])
        return data

    def variables(self):
        return self.varname

    def numRows(self):
        return self.maxrow+1

    def numColumns(self):
        return len(self.header)

    def __len__(self):
        return self.numRows()

    def __getitem__(self, key):
        if type(key) != type(1): return None
        if key < -(self.maxrow+1): return None
        if key > self.maxrow: return None
        return Row(key, self.varmap, self.data[key])
#------------------------------------------------------------------------------
class Buffer:

    def __init__(self, buffer, buffermap, variable):
        self.buffer = buffer
        self.buffermap = buffermap
        self.variable = variable

    def __getattr__(self, variable):
        if self.buffermap.has_key(variable):
            jj = self.buffermap[variable]
            return self.buffer[jj].__getattribute__(variable)
        else:
            raise AttributeError(variable)

    def __call__(self, variable):
        return self.__getattr__(variable)
    
    def __str__(self):
        s = 'Event variables:\n'
        for tname, name, maxcount in self.variable:
            s += "  %-12s %-24s:\t" % (tname, name)
            s += "%s\n" % self.__getattr__(name)
        return s

class Ntuple:
    # "self" is Python's equivalent of the "this" pointer in C++
    # self points to the memory allocated for the object

    def __init__(self, filename, treename, firstrow=0, nrows=None):
        from string import find
        
        # cache inputs
        self.status = 0

        from random import randint
        self.postfix  = randint(1, 1000000)
        if type(filename) == type(""):
            self.filename = [filename]
        else:
            self.filename = filename

        self.currentTreeNumber = -1
        self.treename = treename
        self.nrows = nrows

        # make sure all files exist
        fnames = []
        for fname in self.filename:
            if not os.path.exists(fname):
                print "** Ntuple *** "\
                      "root file %s not found" % fname
                sys.exit(0)

            # check file size
            size = os.path.getsize(fname)
            if size == 0:
                print "== ZERO length == %s" % fname
                continue
            else:
                fnames.append(fname)

        self.filename = fnames
        if len(self.filename) == 0:
            self.status = -3
            return

        # create a chain of files
        self.chain = ROOT.TChain(treename)
        if not self.chain:
            print "*** Ntuple *** can't create chain %s" % treename
            sys.exit(0)

        # ------------------------------------
        # add files to chain
        # ------------------------------------
        for fname in self.filename:
            try:
                self.chain.Add(fname)
            except:
                print "*** problem with file %s" % fname

        self.entries = self.chain.GetEntries()
        self.tree = self.chain
        tree = self.tree

        if tree == None:
            print "** problem accessing Tree - "\
              "perhaps the name %s is wrong?" % treename
            sys.exit(0)

        # get names of variables from root file
        branches  = tree.GetListOfBranches()

        # get number of variables
        try:
            nbranches = branches.GetEntries()
        except:
            print "** ====>  problem accessing branches\n"
            self.status = -1
            return

        bnamemap = {}        
        self.vars = []
        for i in xrange(nbranches):
            # get the ith branch (aka variable)
            bname = branches[i].GetName()
            			
            # just in case, check for duplicates!
            if bnamemap.has_key(bname):
                print '** duplicate branch name: %s' % bname
                continue
            else:
                bnamemap[bname] = 1
            
            # assume one leaf/branch
            # Get all leaves associated with this branch
            leaves = branches[i].GetListOfLeaves()
            if leaves == None:
                print "No leaves found!"
                sys.exit(0)

            leaf = leaves[0]
            if leaf == None:
                print "No leaf found"
                sys.exit(0)

            leafname = leaf.GetName()

            # get leaf type (int, float, double, etc.)
            tname = leaf.GetTypeName()

            if find(tname, 'vector') > -1:
                #print '**skipping %s\t%s for now' % (tname, name)
                continue
            
            # check for leaf counter
            flag = Long(0)
            leafcounter = leaf.GetLeafCounter(flag)
            if leafcounter:
                maxcount = leafcounter.GetMaximum()
            else:
                maxcount = leaf.GetLen()

            # store type and variable name
            self.vars.append( (tname, bname, maxcount) )
                
        nlen = len(self.vars)

        # create a map of variable name to column number
        self.varmap = {}
        for ind, var in enumerate(self.vars):
            self.varmap[var] = ind

        nentries = self.entries
        if self.nrows != None:
            self.entries = min(self.nrows, nentries) 
        else:			
            self.entries = nentries
        # ------------------------------------
        # set up branches as a struct
        # Root has a limit on how long a
        # string it can cope with, so split
        # into multiple strings
        # ------------------------------------

        bufferCount = 0
        newBuffer = True
        rec = ""
        bufferName = ""
        maxlength  = 2000
        self.buffermap  = {}
        self.buffer = []

        self.varname= []
        for count, (tname, name, maxcount) in enumerate(self.vars):
            self.varname.append((name, maxcount))
            
            # keep track of map from variable name to buffer count
            self.buffermap[name] = bufferCount

            if newBuffer:
                newBuffer = False
                bufferName = "S%d_%d" % (self.postfix, bufferCount)
                rec = "struct %s {" % bufferName

            if find(tname, 'vector') > -1:
                # special handling for vectors
                continue
            
            elif maxcount == 1:
                rec += "%s %s;" % (tname, name)
            else:				
                rec += "%s %s[%d];" % (tname, name, maxcount)

            if (len(rec) > maxlength) or \
                   (count >= len(self.vars)-1):
                rec += "};"
                newBuffer = True

                # evaluate it

                ROOT.gROOT.ProcessLine(rec)

                # now import struct
                exec("from ROOT import %s" % bufferName)

                # add to list of buffers
                self.buffer.append(eval("%s()" % bufferName))

                # remember to update buffer count
                bufferCount += 1

        # create a generic event object
        self.event = Buffer(self.buffer, self.buffermap, self.vars)

        # Now that addresses are stable, give address of each variable
        for tname, name, maxcount in self.vars:
            jj = self.buffermap[name]
            tree.SetBranchAddress(name, ROOT.AddressOf(self.buffer[jj], name))

        self.status = 0
        # initialize row number
        self.row = firstrow

    # destructor
    def __del__(self):
        pass

    def close(self):
        pass
    
    def size(self):
        return int(self.entries)

    def numEntries(self):
        return int(self.entries)

    def __len__(self):
        return int(self.entries)
    
    def read(self, row):
        localentry = self.chain.LoadTree(row)
        if self.chain.GetTreeNumber() != self.currentTreeNumber:
            self.currentTreeNumber = self.chain.GetTreeNumber()
            # Update branch addresses
            self.tree  = self.chain.GetTree()
            for tname, name, maxcount in self.vars:
                jj = self.buffermap[name]
                self.tree.SetBranchAddress(name,
                                           ROOT.AddressOf(self.buffer[jj], name))

        self.tree.GetEntry(localentry)

    def treeNumber(self):
        return (self.currentTreeNumber,
                self.filename[self.currentTreeNumber])

    def good(self):
        return self.status == 0

    def get(self, variable):
        if self.buffermap.has_key(variable):
            jj = self.buffermap[variable]
            return self.buffer[jj].__getattribute__(variable)
        else:
            return None

    def __call__(self, variable):
        return self.get(variable)
    
    def __str__(self):
        rec = ''
        for ii, (tname, bname, maxcount) in enumerate(self.vars):
            rec += "\t%4d\t%-12s\t%-32s\t%d\n" % (ii, tname, bname, maxcount)
        return rec

    def ls(self):
        print self.__str__()

    def variables(self):
        return self.varname
    
    # Implement Python iterator protocol
    def __iter__(self):
        return self

    def next(self):
        if self.row > self.entries-1:
            self.row = 0
            raise StopIteration
        else:
            self.read(self.row)
            self.row += 1
            return self.event
#------------------------------------------------------------------------------
class Node:
    def __init__(self,
                 left, right, selector, cutValue,
                 cutType, nodeType, purity, response):
        self.left  = left            # left node
        self.right = right           # right node
        self.selector = selector     # variable index
        self.cutValue = cutValue     # variable value
        self.cutType  = cutType
        # -1, 0, 1 (bkg leaf, internal node, signal leaf)              
        self.nodeType = nodeType     
        self.purity   = purity
        self.response = response

    def __del__(self):
        pass

    # test event if it decends the tree at this node to the right
    def goesRight(self, inputValues):
        result = inputValues[self.selector] > self.cutValue
        if self.cutType:
            return result # the cuts are selecting Signal 
        else:
            return not result

    # test event if it descends the tree at this node to the left 
    def goesLeft (self, inputValues):
        if not self.goesRight(inputValues):
            return True
        else:
            return False
       
    def getRight(self):
        return self.right

    def getLeft(self):
        return self.left   

    def getPurity(self):
        return self.purity
    
    def getNodeType(self):
        return self.nodeType
    
    def getResponse(self):
        return self.response

#-----------------------------------------------------------------------------
class BDT:
    def __init__(self, filename, normweights=True):        
        import re
        from os import path
        from sys import exit
        if not path.exists(filename):
            print '** BDT ** error ** cannot open file %s' % filename
            exit()

        self.normweights = normweights
        
        record = open(filename).read()
        gettree = re.compile('^  [/][/] itree[^)]+[)];[^;]+;', re.M)
        recs = gettree.findall(record)

        # get names of input variables
        getinputvars = re.compile('(?<=const char\* inputVars\[\] = ).*')
        getvars      = re.compile('(?<= ").*?(?=",)|(?<= ").*?(?=" )')
        self.varnames= getvars.findall(getinputvars.findall(record)[0])        
        self.weights = []
        self.forest  = []

        for index, record in enumerate(recs):
            record = replace(record, 'NN(', 'Node(')
            record = replace(record,'  //','#')
            record = replace(record,'  fBoostWeights.push_back',
                                 'self.weights.append')
            record = replace(record,'  fForest.push_back',
                                 'self.forest.append')
            record = replace(record,';','')
            exec(record)
                
    def __del__(self):
        pass

    def __call__(self, inputValues, numTrees=-1):

        totalTrees = len(self.forest)
        if numTrees > 0:
            ntrees = min(numTrees, totalTrees)
        else:
            ntrees = totalTrees

        value = 0.0
        norm  = 0.0
        for itree in xrange(ntrees):
            current = self.forest[itree]
            while current.getNodeType() == 0:
                if current.goesRight(inputValues):
                    current = current.getRight()
                else:
                    current = current.getLeft()
            value += self.weights[itree] * current.getNodeType()
            norm  += self.weights[itree]
        if self.normweights:
            value /= norm
        else:
            value = 1.0/(1 + exp(-2*value))
        return value

    def __len__(self):
        return len(self.forest)
    
    def variables(self):
        return self.varnames
    
    def printTree(self, itree, depth=0, which=0, node=None):
        if which == 0:
            node = self.forest[itree]
            print "tree number %d\tweight = %10.3e" % \
              (itree, self.weights[itree])

        if depth > 100:  return
        if node == 0:    return
        if node == None: return
            
        selector = node.selector
        if selector > -1:
            name = self.varnames[node.selector]
        else:
            name = 'LEAF '

        value = node.cutValue
        if which == 0:
            nodedir = 'root '
        elif which < 0:
            nodedir = 'left '
        else:
            nodedir = 'right'

        if node.getNodeType() < 0:
            nodetype = 'B'
            value = node.getPurity()
        elif node.getNodeType() > 0:
            nodetype = 'S'
            value = node.getPurity()

        depthstr = "  "*(depth+1)
        print "%s %10s %10s\t%10.2f" % (depthstr, which, name, value)
        
        depth += 1
        self.printTree(itree, depth, -1, node.left)
        self.printTree(itree, depth,  1, node.right)

        
    def ranking(self, ntrees=-1):
        self.countname = {}
        if ntrees <= 0:
            maxtrees = len(self.forest)
        else:
            maxtrees = min(ntrees, len(self.forest))
        for itree in xrange(maxtrees):
            self.__rank(itree)

        recs = []
        total = 0.0
        for (key, value) in self.countname.items():
            total += value
            recs.append([value, key])
        for index, (value, key) in enumerate(recs):
            recs[index] = (recs[index][0]/total, recs[index][1])
        recs.sort()
        recs.reverse()
        return recs
            
    def __rank(self, itree, depth=0, which=0, node=None):
        if which == 0:
            node = self.forest[itree]
            
        if (node == 0) or (node == None) or (node.selector < 0):
            return
        
        name = self.varnames[node.selector]
        if not self.countname.has_key(name):
            self.countname[name] = 0
        self.countname[name] += 1.0
        depth += 1
        self.__rank(itree, depth, -1, node.left)
        self.__rank(itree, depth,  1, node.right)

        
    def weight(self, itree):
        if itree >=0 and itree < len(self.weights):
            return self.weights[itree]
        else:
            return -1
        
    def summedWeights(self):
        return sum(self.weights)
    
    def setValueType(self, which=0):
        self.valueType = which
        
    # plot only valid for 2-D applications
    def plot2d(self, itree, hname, xtitle, ytitle,
             xmin, xmax, ymin, ymax, useValue=False,
             node=None):
        
        if node == None:
            hname = "%s%5.5d" % (hname, itree)
            self.hplot = ROOT.TH2Poly(hname, "", xmin, xmax, ymin, ymax)
            self.hplot.GetXaxis().SetTitle(xtitle)
            self.hplot.GetYaxis().SetTitle(ytitle)
            self.hplot.GetYaxis().SetTitleOffset(1.6)
            self.hplot.SetNdivisions(505, "X")
            self.hplot.SetNdivisions(505, "Y")
            self.binNumber = 0
            node = self.forest[itree]

        if node == None:
            print "*** node is None - shouldn't happen ***"
            sys.exit(0)

        if self.hplot == None:
            print "*** hplot is None - shouldn't happen ***"
            sys.exit(0)

        if node == 0:
            return self.hplot
 
        self.binNumber += 1
        if self.binNumber > 20:
            print "*** lost in trees ***"
            sys.exit(0)
            
        if useValue:
            weight = node.nodeType        
        else:
            weight =-self.binNumber
            
        self.hplot.AddBin(xmin, ymin, xmax, ymax)            
        self.hplot.SetBinContent(self.binNumber, weight)

        if node.selector < 0:
            return self.hplot
        
        value = node.cutValue
        if node.selector == 0:
            # left
            xmax1= xmax
            xmax = value
            self.plot2d(itree, hname, xtitle, ytitle,
                        xmin, xmax, ymin, ymax, useValue, node.left)
            # right
            xmax = xmax1
            xmin = value
            self.plot2d(itree, hname, xtitle, ytitle,
                        xmin, xmax, ymin, ymax, useValue, node.right)
        else:
            # left
            ymax1 = ymax
            ymax = value
            self.plot2d(itree, hname, xtitle, ytitle,
                        xmin, xmax, ymin, ymax, useValue, node.left)
            # right
            ymax = ymax1
            ymin = value
            self.plot2d(itree, hname, xtitle, ytitle,
                        xmin, xmax, ymin, ymax, useValue, node.right)
        return self.hplot

