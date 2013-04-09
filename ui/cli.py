# -*-coding:utf-8-*-
import sys, ConfigParser
from collections import OrderedDict

def IRCLI():
    """
    image render command line interface
    """
    algos = ['CDIR']
    menu = "image rendering algorithms:\n" 
    for i in range(len(algos)):
        menu += "%d. %s " %(i+1, algos[i])
    menu += "\n"
    
    prompt = [
        ("algorithm", menu + "please choose which one to use", 1, len(algos)),
        ("minstrokesize", "minimum stroke size", 2, 32),
        ("maxstrokesize", "maximum stroke size", 2, 32),
        ("strokesizestep", "stroke size step", 1, 30),
        ("minstrokelen", "minimum stroke length", 1, 32),
        ("maxstrokelen", "maximum stroke length", 1, 32),
        ]
    params = OrderedDict()

    for item in prompt:
        key = item[0]
        msg = item[1]
        minVal = item[2]
        maxVal = item[3]
        choice = -1
        while choice < minVal or choice > maxVal:
            choice = input("%s (%d ~ %d): " %(msg, minVal, maxVal))
        params[key] = choice

    params["algorithm"] = algos[params["algorithm"]-1]
    return params

def VRCLI():
    """
    video render command line interface
    """
    algos = ['LAOF']
    menu = "video rendering algorithms:\n" 
    for i in range(len(algos)):
        menu += "%d. %s " %(i+1, algos[i])
    menu += "\n"
    
    prompt = [
        ("algorithm", menu + "please choose which one to use", 1, len(algos)),
        ("spatial", "spatial smooth factor"),
        ("temporal", "temporal smooth factor"),
        ("segment", "segment  similarity"),
        ]
    params = OrderedDict()

    for item in prompt:
        key = item[0]
        msg = item[1]

        if len(item) > 2:
            minVal = item[2]
            maxVal = item[3]
            choice = -1
            while choice < 0 or choice > maxVal:
                choice = input("%s (%d ~ %d): " %(msg, minVal, maxVal))
        else:
            choice = input("%s (> 0): " %msg)
            
        params[key] = choice

    params["algorithm"] = algos[params["algorithm"]-1]
    return params

def dict2conf(params):
    conf = ConfigParser.ConfigParser()
    
    for sec, options in params.items():
        conf.add_section(sec)
        for option, value in options.items():
            conf.set(sec, option, value)
    return conf

def conf2dict(conf):
    params = OrderedDict()
    for sec in conf.sections():
        params[sec] = {}
        for option, value in conf.items(sec):
            params[sec][option] = value
    return params

if __name__ == "__main__":
    if len(sys.argv) == 1:
        infile = raw_input("please input video to process: ")
        outfile = raw_input("please input path to save: ")
        ir = IRCLI()
        vr = VRCLI()
        params = OrderedDict([
            ("Files", OrderedDict([("infile",infile), ("outfile",outfile)])),
            ("ImageRender", ir), ("VideoRender", vr)])
        config = dict2conf(params)
        with open("params.conf", "w") as fp:
            config.write(fp)
    else:
        config = ConfigParser.ConfigParser()
        print "reading config file %s... done" %sys.argv[1]
        with open(sys.argv[1]) as fp:
            config.readfp(fp)
        params = conf2dict(config)
        infile = params["Files"]["infile"]
        outfile = params["Files"]["outfile"]

    print "converting %s into frames... done" %infile
    print "calculating dense correspondences... done"
    print "stroke propagation... done"
    print "saving to %s... done" %outfile
