#!/usr/bin/python
#
# Copyright 2012 Cyril Plisko. All rights reserved.
# Use is subject to license terms.
#

import argparse
import cStringIO
import os.path
import time
import subprocess as sp

VERSION = '1'

MODPARAMPREF = '/sys/module/zfs/parameters'
ZFSMODPARAMS = ['/sbin/modinfo', '-F', 'parm', 'zfs']
ARCSTATSPATH = '/proc/spl/kstat/zfs/arcstats'
ZILSTATSPATH = '/proc/spl/kstat/zfs/zil'

ARCSUMMARY = '''
ARC Size
--------
\tMinimal\t\t{c_min!s}
\tCurrent\t\t{size!s}
\tTarget\t\t{c!s}
\tMaximal\t\t{c_max!s}

ARC Size Breakdown
------------------
\tMRU Size\t{p!s}
\tMFU Size\t{mfu!s}

ARC Efficiency
--------------
\tHit ratio\t{ratio:.2%}

ARC Meta
--------

\tARC Metadata\t{arc_meta_used!s}

'''

L2ARCSUMMARY = '''
L2ARC
-----
\tSize\t\t{l2_size!s}
\tHeader Size\t{l2_hdr_size!s} (ARC overhead)
\tHit ratio\t{l2_ratio:.2%}
'''

METRIC_NAMES = {'hps': 'hit/s', 'mps': 'miss/s',
    'ddhps': 'datah/s', 'ddmps': 'datam/s',
    'dmdhps': 'metah/s', 'dmdmps': 'metam/s',
    'pdhps': 'datah/s', 'pdmps': 'datam/s',
    'pmdhps': 'metah/s', 'pmdmps': 'metam/s',
    'c': 'deltac', 'p': 'deltap',
    }

OUT_HEADER = '{total:^16}{demand:^32}{prefetch:^32}'
OUT_TOTAL = '{hps:>8}{mps:>8}'
OUT_DEMAND = '{ddhps:>8}{ddmps:>8}{dmdhps:>8}{dmdmps:>8}'
OUT_PREFETCH = '{pdhps:>8}{pdmps:>8}{pmdhps:>8}{pmdmps:>8}'
OUT_ARC = '{c!s:>8}{p!s:>8}'
OUT_FORMAT = OUT_TOTAL + OUT_DEMAND + OUT_PREFETCH + OUT_ARC


def humanize(number):
    number = int(number)
    if abs(number) < 1024:
        return '{n} B'.format(n=number)
    elif abs(number) < 1024 * 1024:
        return '{n} KiB'.format(n=number / 1024)
    elif abs(number) < 1024 * 1024 * 1024:
        return '{n} MiB'.format(n=number / (1024 * 1024))
    elif abs(number) < 1024 * 1024 * 1024 * 1024:
        return '{n} GiB'.format(n=number / ( 1024 * 1024 * 1024))
    elif abs(number) < 1024 * 1024 * 1024 * 1024 * 1024:
        return '{n} TiB'.format(n=number / ( 1024 * 1024 * 1024 * 1024))
    else:
        return '{n} PiB'.format(n=number / ( 1024 * 1024 * 1024 * 1024 * 1024))


def humanize_dict(d):
    res = dict()
    for i in d:
        res[i] = humanize(d[i])
    return res


def get_arcstats():
    return file(ARCSTATSPATH).read()


def zfsparams():
    cmd = sp.Popen(ZFSMODPARAMS, stdout=sp.PIPE)
    paramlist = cmd.communicate()[0].splitlines()
    params = dict([i.split(':') for i in paramlist])
    out = 'ZFS Module Parameters\n---------------------\n'
    for parm in params:
        val = open(os.path.join(MODPARAMPREF, parm)).read().strip()
        out += '{name:40}{value:>20}  ({parm})\n'.format(name=params[parm], parm=parm, value=val)
    return out


class Integer(int):
    def __str__(self):
        return humanize(self.numerator)


class kstat():
    KSTATBASE = '/proc/spl/kstat'
    # kstat types
    KSTAT_TYPE_RAW = 0
    KSTAT_TYPE_NAMED = 1
    KSTAT_TYPE_INTR = 2
    KSTAT_TYPE_IO = 3
    KSTAT_TYPE_TIMER = 4
    # kstat data types
    KSTAT_DATA_CHAR = 0
    KSTAT_DATA_UINT64 = 4

    def __init__(self, module, name):
        self._file = os.path.join(self.KSTATBASE, module, name)
        self._kstat = dict()
        lines = open(self._file).readlines()
        (self._kid, self._type, self._flags, self._ndata, self._data_size,
            self._crtime, self._snaptime) = [int(i, 0) for i in lines[0].split()]
        if self._type == self.KSTAT_TYPE_RAW:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_TYPE_NAMED:
            self._init_named(lines)
        else:
            raise NotImplementedError(self._type)

    def _init_named(self, lines):
        for line in lines[2:]:
            name, dtype, value = line.split()
            if dtype == self.KSTAT_DATA_CHAR:
                pass
            if dtype == self.KSTAT_DATA_UINT64:
                self._kstat[kname] = Integer(kvalue)
            else:
                pass

    def __str__(self):
        return str(self.__dict__)


class arcstats():
    def __init__(self, text=None):
        self._kstat = dict()
        self._arcstats = dict()
        if text:
            lines = text.splitlines()
            n0, n1, n2, n3, n4, n5, tstamp = lines[0].split()
            self._kstat['tstamp'] = int(tstamp)
            for line in lines[2:]:
                name, unused, value = line.split()
                self._arcstats[name] = Integer(value)
                #print 'Processed', name, value

    def __sub__(self, other):
        diff = arcstats()
        for item in self._kstat:
            diff._kstat[item] = self._kstat[item] - other._kstat[item]
        for item in self._arcstats:
            diff._arcstats[item] = Integer(self._arcstats[item] - other._arcstats[item])
        return diff

    def __str__(self):
        return str(self._arcstats)

    def __repr__(self):
        return str(self._arcstats)

    def debug(self):
        return self._kstat

    def summary(self):
        stats = self._arcstats.copy()
        hits = float(self._arcstats['hits'])
        misses = float(self._arcstats['misses'])
        l2_hits = float(self._arcstats['l2_hits'])
        l2_misses = float(self._arcstats['l2_misses'])
        stats['ratio'] = hits / (hits + misses)
        stats['mfu'] = Integer(stats['c'] - stats['p'])
        result = ARCSUMMARY.format(**stats)
        if self._arcstats.get('l2_size'):
            stats['l2_ratio'] = l2_hits / (l2_hits + l2_misses)
            result += L2ARCSUMMARY.format(**stats)
        return result

    def compute(self):
        delta = self._kstat['tstamp'] / 1000000000
        raw = self._arcstats.copy()
        raw['hps'] = self._arcstats['hits'] / delta
        raw['mps'] = self._arcstats['misses'] / delta
        raw['ddhps'] = self._arcstats['demand_data_hits'] / delta
        raw['ddmps'] = self._arcstats['demand_data_misses'] / delta
        raw['dmdhps'] = self._arcstats['demand_metadata_hits'] / delta
        raw['dmdmps'] = self._arcstats['demand_metadata_misses'] / delta
        raw['pdhps'] = self._arcstats['prefetch_data_hits'] / delta
        raw['pdmps'] = self._arcstats['prefetch_data_misses'] / delta
        raw['pmdhps'] = self._arcstats['prefetch_metadata_hits'] / delta
        raw['pmdmps'] = self._arcstats['prefetch_metadata_misses'] / delta
        return raw


class zil(kstat):
    def __init__(self):
        kstat.__init__(self, 'zfs', 'zil')

    def summary(self):
        return self.summary()


def headers():
    header = '\n'
    groups = dict(total='  TOTAL', demand='DEMAND', prefetch='PREFETCH')
    header += OUT_HEADER.format(**groups)
    header += '\n'
    header += OUT_FORMAT.format(**METRIC_NAMES)
    print header


def data(obj):
    raw = obj.compute()
    #print raw
    print OUT_FORMAT.format(**raw)


def cycle(interval, count):
    step = 1 if count else 0
    count = count if count else 1
    cur = arcstats(get_arcstats())
    time.sleep(1)
    lines = 0
    while count:
        count -= step
        prev = cur
        cur = arcstats(get_arcstats())
        delta = cur - prev
        if lines % 20 == 0:
            headers()
        lines += 1
        data(delta)
        if count:
            time.sleep(interval)


def execute(args):
    if args.debug:
        print arcstats(get_arcstats()).debug()
        as1 = arcstats(get_arcstats())
        time.sleep(1)
        as2 = arcstats(get_arcstats())
        print (as2 - as1)
    elif args.summary:
        print arcstats(get_arcstats()).summary()
    elif args.parm:
        print zfsparams()
    elif args.zil:
        print zil()
    else:
        cycle(args.interval, args.count)

    #pp.pprint(statsnap)


def main():
    parser = argparse.ArgumentParser(
        description='Report various ZFS ARC statistics',
        version=VERSION)
    parser.add_argument('-d', '--debug',
        help='help debug this tool', action='store_true')
    parser.add_argument('-p', '--parm',
        help='print ZFS module parameters', action='store_true')
    parser.add_argument('-s', '--summary',
        help='print ARC summary information', action='store_true')
    parser.add_argument('-z', '--zil',
        help='print ZIL statistics', action='store_true')
    parser.add_argument('interval',
        help='seconds between probes', type=int, nargs='?', default=1)
    parser.add_argument('count',
        help='number of probes to run', type=int, nargs='?', default=0)

    args = parser.parse_args()
    #print args
    execute(args)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
