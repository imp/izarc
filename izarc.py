#!/usr/bin/python
#
# Copyright 2012 Cyril Plisko. All rights reserved.
# Use is subject to license terms.
#

import argparse
import os.path
import time
import pprint as pp
import subprocess as sp

VERSION = '3'

NANOSEC = 1000000000
MODPARAMPREF = '/sys/module/zfs/parameters'
ZFSMODPARAMS = ['/sbin/modinfo', '-F', 'parm', 'zfs']

ARCSUMMARY = '''
ARC Size
--------
\tMinimal\t\t{c_min!s}
\tCurrent\t\t{size!s}
\tTarget\t\t{c!s}
\tMaximal\t\t{c_max!s}

ARC Size Breakdown
------------------
\tMRU Size               \t\t{p!s}
\tMFU Size               \t\t{mfu!s}
\tAnonymous Size         \t\t{anon_size!s}
\tCurrent MRU Size       \t\t{mru_size!s}
\tCurrent Ghost MRU Size\t\t{mru_ghost_size!s}
\tCurrent MFU Size       \t\t{mfu_size!s}
\tCurrent Ghost MFU Size\t\t{mfu_ghost_size!s}

ARC Efficiency
--------------
\tHit ratio\t{ratio:.2%}

ARC Hits Breakdown
------------------
\tMRU Hits      \t\t{mru_hits_ratio:.2%}
\tMFU Hits      \t\t{mfu_hits_ratio:.2%}
\tMRU Ghost Hits\t\t{mru_ghost_hits_ratio:.2%}
\tMFU Ghost Hits\t\t{mfu_ghost_hits_ratio:.2%}

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
    'cps': 'commit/s', 'wps': 'wrcomm/s',
    'itxs': 'itx/s',
    'imnbs': 'imsnb/s', 'imncs': 'imsnc/s',
    'imsbs': 'imssb/s', 'imscs': 'imssc/s',
    'incbs': 'incb/s', 'inccs': 'incc/s',
    'icbs': 'icb/s', 'iccs': 'icc/s',
    'anons': 'dltsize', 'adevc': 'devict/s', 'amevc': 'mevict/s',
    'mrus': 'dltsize', 'mruh': 'hits/s', 'mrudevc': 'devict/s', 'mrumevc': 'mevict/s',
    'mfus': 'dltsize', 'mfuh': 'hits/s', 'mfudevc': 'devict/s', 'mfumevc': 'mevict/s',
    'gmrus': 'dltsize', 'gmruh': 'hits/s', 'gmrudevc': 'devict/s', 'gmrumevc': 'mevict/s',
    'gmfus': 'dltsize', 'gmfuh': 'hits/s', 'gmfudevc': 'devict/s', 'gmfumevc': 'mevict/s',
    'l2hps': 'hit/s', 'l2mps': 'miss/s',
    'l2size': 'dltsize', 'l2hdrsize': 'dlthdrsize',
    'l2read': 'read/s', 'l2write': 'write/s', 'l2sent': 'sent/s',
    'l2whm': 'writehdrm/s', 'l2fow': 'freeonwrite/s',
    }

HEADER_NAMES = {'total': '  TOTAL', 'demand': 'DEMAND', 'prefetch': 'PREFETCH',
    'arc': 'ARC SIZE', 'transactions': 'TRANSACTIONS', 'copy': 'DATA COPY', 'anon': 'ANON',
    'mru': 'MRU', 'mfu': 'MFU', 'gmru': 'GHOST MRU', 'gmfu': 'GHOST MFU'}

ARC_HEADER = '{total:^16}{demand:^32}{prefetch:^32}{arc:^16}'
OUT_TOTAL = '{hps:>8}{mps:>8}'
OUT_DEMAND = '{ddhps:>8}{ddmps:>8}{dmdhps:>8}{dmdmps:>8}'
OUT_PREFETCH = '{pdhps:>8}{pdmps:>8}{pmdhps:>8}{pmdmps:>8}'
OUT_ARC = '{c!s:>8}{p!s:>8}'
ARC_FORMAT = OUT_TOTAL + OUT_DEMAND + OUT_PREFETCH + OUT_ARC

EXTENDEDARC_HEADER = '{anon:^28}{mru:^36}{mfu:^36}{gmru:^36}{gmfu:^36}'
OUT_ANON = '{anons!s:>8}{adevc:>10}{amevc:>10}'
OUT_MRU = '{mrus!s:>8}{mruh:>8}{mrudevc:>10}{mrumevc:>10}'
OUT_MFU = '{mfus!s:>8}{mfuh:>8}{mfudevc:>10}{mfumevc:>10}'
OUT_GMRU = '{gmrus!s:>8}{gmruh:>8}{gmrudevc:>10}{gmrumevc:>10}'
OUT_GMFU = '{gmfus!s:>8}{gmfuh:>8}{gmfudevc:>10}{gmfumevc:>10}'
EXTENDEDARC_FORMAT = OUT_ANON + OUT_MRU + OUT_MFU + OUT_GMRU + OUT_GMFU

L2ARC_SIZE = '{l2size!s:>14}{l2hdrsize!s:>14}'
L2ARC_IO = '{l2hps:>8}{l2mps:>8}{l2read:>8}{l2write:>8}{l2sent:>8}'
L2ARC_MISC = '{l2whm:>14}{l2fow:>14}'
L2ARC_FORMAT = L2ARC_SIZE + L2ARC_IO + L2ARC_MISC

ZIL_HEADER = '{total:^20}{transactions:^40}{copy:^32}'
ZIL_COMMITS = '{cps:>10}{wps:>10}'
ZIL_ITX = '{itxs:>8}{imnbs!s:>8}{imncs:>8}{imsbs!s:>8}{imscs:>8}'
ZIL_COPIES = '{incbs!s:>8}{inccs:>8}{icbs!s:>8}{iccs:>8}'
ZIL_FORMAT = ZIL_COMMITS + ZIL_ITX + ZIL_COPIES


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


def zfsparams():
    cmd = sp.Popen(ZFSMODPARAMS, stdout=sp.PIPE)
    paramlist = cmd.communicate()[0].splitlines()
    params = dict([i.split(':') for i in paramlist])
    out = 'ZFS Module Parameters\n---------------------\n'
    for parm in params:
        val = open(os.path.join(MODPARAMPREF, parm)).read().strip()
        out += '{name:40}{value:>20}  ({parm})\n'.format(name=params[parm], parm=parm, value=val)
    return out


class Integer(long):
    def __str__(self):
        return humanize(self.numerator)


class kstat(object):
    KSTATBASE = '/proc/spl/kstat'
    # kstat types
    KSTAT_TYPE_RAW = 0
    KSTAT_TYPE_NAMED = 1
    KSTAT_TYPE_INTR = 2
    KSTAT_TYPE_IO = 3
    KSTAT_TYPE_TIMER = 4
    KSTAT_TYPE_TXG = 5

    # kstat data types
    KSTAT_DATA_CHAR = 0
    KSTAT_DATA_UINT64 = 4

    def __init__(self, module=None, name=None):
        self._kstat = dict()
        if module and name:
            self._file = os.path.join(self.KSTATBASE, module, name)
            self._init_kstat()
        else:
            self._file = ''

    def _init_kstat(self):
        lines = open(self._file).readlines()
        # Parse common kstat data
        (self._kid, self._type, self._flags, self._ndata, self._data_size,
            self._crtime, self._snaptime) = [int(i, 0) for i in lines[0].split()]
        # Parse the rest of kstat data
        if self._type == self.KSTAT_TYPE_RAW:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_TYPE_NAMED:
            self._init_named(lines)
        elif self._type == self.KSTAT_TYPE_INTR:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_TYPE_IO:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_TYPE_TIMER:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_TYPE_TXG:
            raise NotImplementedError(self._type)
        else:
            raise NotImplementedError(self._type)

    def _init_named(self, lines):
        for line in lines[2:]:
            name, data, value = line.split()
            data = int(data)
            if data == self.KSTAT_DATA_CHAR:
                pass
            elif data == self.KSTAT_DATA_UINT64:
                self._kstat[name] = int(value, 0)
            else:
                pass

    def __sub__(self, other):
        diff = self.__class__()
        diff._snaptime = self._snaptime - other._snaptime
        for item in self._kstat:
            diff._kstat[item] = Integer(self._kstat[item] - other._kstat[item])
        return diff

    def __str__(self):
        text = 'kid {self._kid}\n'.format(self=self)
        text = pp.pformat(self._kstat)
        return text


class arcstats(kstat):
    def __init__(self):
        super(arcstats, self).__init__('zfs', 'arcstats')
        self._arcstats = dict()
        for name in self._kstat:
            self._arcstats[name] = Integer(self._kstat[name])

        if False:
            lines = text.splitlines()
            n0, n1, n2, n3, n4, n5, tstamp = lines[0].split()
            self._kstat['tstamp'] = int(tstamp)
            for line in lines[2:]:
                name, unused, value = line.split()
                self._arcstats[name] = Integer(value)

    def __sub__(self, other):
        diff = self.__class__()
        diff._snaptime = self._snaptime - other._snaptime
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
        stats['mru_hits_ratio'] = float(self._arcstats['mru_hits']) / hits
        stats['mfu_hits_ratio'] = float(self._arcstats['mfu_hits']) / hits
        stats['mru_ghost_hits_ratio'] = float(self._arcstats['mru_ghost_hits']) / hits
        stats['mfu_ghost_hits_ratio'] = float(self._arcstats['mfu_ghost_hits']) / hits
        result = ARCSUMMARY.format(**stats)
        if self._arcstats.get('l2_size'):
            stats['l2_ratio'] = l2_hits / (l2_hits + l2_misses)
            result += L2ARCSUMMARY.format(**stats)
        return result

    def compute(self):
        delta = self._snaptime / NANOSEC
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

    def headers(self):
        header = '\n'
        header += ARC_HEADER.format(**HEADER_NAMES)
        header += '\n'
        header += ARC_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return ARC_FORMAT.format(**raw)


class extendarc(arcstats):
    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._arcstats.copy()
        raw['anons'] = self._arcstats['anon_size']
        raw['adevc'] = self._arcstats['anon_evict_data'] / delta
        raw['amevc'] = self._arcstats['anon_evict_metadata'] / delta
        raw['mrus'] = self._arcstats['mru_size']
        raw['mruh'] = self._arcstats['mru_hits'] / delta
        raw['mrudevc'] = self._arcstats['mru_evict_data'] / delta
        raw['mrumevc'] = self._arcstats['mru_evict_metadata'] / delta
        raw['mfus'] = self._arcstats['mfu_size']
        raw['mfuh'] = self._arcstats['mfu_hits'] / delta
        raw['mfudevc'] = self._arcstats['mfu_evict_data'] / delta
        raw['mfumevc'] = self._arcstats['mfu_evict_metadata'] / delta
        raw['gmrus'] = self._arcstats['mru_ghost_size']
        raw['gmruh'] = self._arcstats['mru_ghost_hits'] / delta
        raw['gmrudevc'] = self._arcstats['mru_ghost_evict_data'] / delta
        raw['gmrumevc'] = self._arcstats['mru_ghost_evict_metadata'] / delta
        raw['gmfus'] = self._arcstats['mfu_ghost_size']
        raw['gmfuh'] = self._arcstats['mfu_ghost_hits'] / delta
        raw['gmfudevc'] = self._arcstats['mfu_ghost_evict_data'] / delta
        raw['gmfumevc'] = self._arcstats['mfu_ghost_evict_metadata'] / delta
        return raw

    def headers(self):
        header = '\n'
        header += EXTENDEDARC_HEADER.format(**HEADER_NAMES)
        header += '\n'
        header += EXTENDEDARC_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return EXTENDEDARC_FORMAT.format(**raw)


class l2arc(arcstats):
    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._arcstats.copy()
        raw['l2size'] = self._arcstats['l2_size']
        # L2ARC headers (ARC overhead)
        raw['l2hdrsize'] = self._arcstats['l2_hdr_size']
        raw['l2hps'] = self._arcstats['l2_hits'] / delta
        raw['l2mps'] = self._arcstats['l2_misses'] / delta
        raw['l2read'] = self._arcstats['l2_read_bytes'] / delta
        raw['l2write'] = self._arcstats['l2_write_bytes'] / delta
        raw['l2sent'] = self._arcstats['l2_writes_sent'] / delta
        # This buffer misses out.  It may be in a stage
        # of eviction.  Its ARC_L2_WRITING flag will be
        # left set, denying reads to this buffer.
        raw['l2whm'] = self._arcstats['l2_writes_hdr_miss'] / delta
        # Free the arc data buffer.  If it is an l2arc write in progress,
        # the buffer is placed on l2arc_free_on_write to be freed later.
        raw['l2fow'] = self._arcstats['l2_free_on_write'] / delta
        return raw

    def headers(self):
        header = '\n'
        header += L2ARC_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return L2ARC_FORMAT.format(**raw)


class zil(kstat):
    def __init__(self):
        super(zil, self).__init__('zfs', 'zil')

    def summary(self):
        return self.summary()

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._kstat.copy()
        raw['cps'] = self._kstat['zil_commit_count'] / delta
        raw['wps'] = self._kstat['zil_commit_writer_count'] / delta
        raw['itxs'] = self._kstat['zil_itx_count'] / delta
        raw['ibs'] = Integer(self._kstat['zil_itx_indirect_bytes'] / delta)
        raw['ics'] = self._kstat['zil_itx_indirect_count'] / delta
        raw['imnbs'] = Integer(self._kstat['zil_itx_metaslab_normal_bytes'] / delta)
        raw['imncs'] = self._kstat['zil_itx_metaslab_normal_count'] / delta
        raw['imsbs'] = Integer(self._kstat['zil_itx_metaslab_slog_bytes'] / delta)
        raw['imscs'] = self._kstat['zil_itx_metaslab_slog_count'] / delta
        raw['incbs'] = Integer(self._kstat['zil_itx_needcopy_bytes'] / delta)
        raw['inccs'] = self._kstat['zil_itx_needcopy_count'] / delta
        raw['icbs'] = Integer(self._kstat['zil_itx_copied_bytes'] / delta)
        raw['iccs'] = self._kstat['zil_itx_copied_count'] / delta
        return raw

    def headers(self):
        header = '\n'
        header += ZIL_HEADER.format(**HEADER_NAMES)
        header += '\n'
        header += ZIL_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return ZIL_FORMAT.format(**raw)


def cycle(stats, interval, count):
    step = 1 if count else 0
    count = count if count else 1
    cur = stats()
    time.sleep(1)
    lines = 0
    while count:
        count -= step
        prev = cur
        cur = stats()
        delta = cur - prev
        if lines % 20 == 0:
            print delta.headers()
        lines += 1
        print delta.data()
        if count:
            time.sleep(interval)


def execute(args):
    if args.debug:
        print arcstats().debug()
        as1 = arcstats()
        time.sleep(1)
        as2 = arcstats()
        print (as2 - as1)
    elif args.summary:
        print arcstats().summary()
    elif args.parm:
        print zfsparams()
    elif args.zil:
        if args.debug:
            print zil()
        cycle(zil, args.interval, args.count)
    elif args.l2arc:
        cycle(l2arc, args.interval, args.count)
    elif args.extend:
        cycle(extendarc, args.interval, args.count)
    else:
        cycle(arcstats, args.interval, args.count)


def main():
    parser = argparse.ArgumentParser(
        description='Report various ZFS ARC statistics')
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    parser.add_argument('-d', '--debug',
        help='help debug this tool', action='store_true')
    parser.add_argument('-p', '--parm',
        help='print ZFS module parameters', action='store_true')
    parser.add_argument('-s', '--summary',
        help='print ARC summary information', action='store_true')
    parser.add_argument('-z', '--zil',
        help='print ZIL statistics', action='store_true')
    parser.add_argument('-x', '--extend',
        help='print extended ARC statistics', action='store_true')
    parser.add_argument('-l', '--l2arc',
        help='print L2ARC statistics', action='store_true')
    parser.add_argument('interval',
        help='seconds between probes', type=int, nargs='?', default=1)
    parser.add_argument('count',
        help='number of probes to run', type=int, nargs='?', default=0)

    args = parser.parse_args()
    execute(args)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
