#!/usr/bin/python
#
# Copyright 2012 Cyril Plisko. All rights reserved.
# Use is subject to license terms.
#

import sys
import argparse
import os.path
import time
import pprint as pp
import subprocess as sp
from copy import deepcopy
from decimal import Decimal, getcontext, InvalidOperation, DivisionByZero
from datetime import datetime
from glob import glob

VERSION = '3'

NANOSEC = 1000000000
MODPARAMPREF = '/sys/module/zfs/parameters'
ZFSMODPARAMS = ['/sbin/modinfo', '-F', 'parm', 'zfs']
ZFSVERSION = ['rpm', '-q', 'zfs']

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
    # zil metrics
    'cps': 'commit/s', 'wps': 'wrcomm/s',
    'itxs': 'itx/s',
    'imnbs': 'imsnb/s', 'imncs': 'imsnc/s',
    'imsbs': 'imssb/s', 'imscs': 'imssc/s',
    'incbs': 'incb/s', 'inccs': 'incc/s',
    'iibs': 'iib/s', 'iics': 'iic/s',
    'icbs': 'icb/s', 'iccs': 'icc/s',
    # extend arc metrics
    'anons': 'dltsize', 'adevc': 'devict', 'amevc': 'mevict',
    'mrus': 'dltsize', 'mruh': 'hits/s', 'mrudevc': 'devict', 'mrumevc': 'mevict',
    'mfus': 'dltsize', 'mfuh': 'hits/s', 'mfudevc': 'devict', 'mfumevc': 'mevict',
    'gmrus': 'dltsize', 'gmruh': 'hits/s', 'gmrudevc': 'devict', 'gmrumevc': 'mevict',
    'gmfus': 'dltsize', 'gmfuh': 'hits/s', 'gmfudevc': 'devict', 'gmfumevc': 'mevict',
    # l2arc metrics
    'l2hps': 'hit/s', 'l2mps': 'miss/s',
    'l2size': 'dltsize', 'l2hdrsize': 'dlthdrsize',
    'l2read': 'read/s', 'l2write': 'write/s',
    # prefetch metrics
    'pfhps': 'hit/s', 'pfmps': 'miss/s',
    'pfchps': 'chit/s', 'pfcmps': 'cmiss/s',
    'pfshps': 'shit/s', 'pfsmps': 'smiss/s',
    'pfrs': 'rsuc/s', 'pfrf': 'rfail/s', 'pfbs': 'sbogus/s',
    'pfsr': 'sreset/s', 'pfsnr': 'snoreset/s',
    # io metrics
    'nreads': 'read/s', 'nwrittens': 'write/s', 'readss': 'read/s', 'writess': 'write/s',
    'awaitq': 'awaitq', 'waitpct': 'waitq%',
    'arunq': 'arunq', 'runpct': 'runq%',
    # txgs metrics
    'txg': 'txg', 'birth': 'birth', 'state': 'state',
    'nreserved': 'nreserved', 'nread': 'nread', 'nwritten': 'nwritten', 'reads': 'reads', 'writes': 'writes',
    'otime': 'otime', 'qtime': 'qtime', 'wtime': 'wtime', 'stime': 'stime',
    # tx_assign metrics
    'atime': 'atime/s',
    }


HEADER_NAMES = {'total': '  TOTAL', 'demand': 'DEMAND', 'prefetch': 'PREFETCH',
    'arc': 'ARC SIZE', 'transactions': 'TRANSACTIONS', 'copy': 'DATA COPY', 'anon': 'ANON',
    'mru': 'MRU', 'mfu': 'MFU', 'gmru': 'GHOST MRU', 'gmfu': 'GHOST MFU', 'bandwidth': 'BANDWIDTH',
    'operations': 'OPERATIONS', 'queue': 'QUEUE'}

ARC_HEADER = '{total:^16}{demand:^32}{prefetch:^32}{arc:^20}'
OUT_TOTAL = '{hps:>8}{mps:>8}'
OUT_DEMAND = '{ddhps:>8}{ddmps:>8}{dmdhps:>8}{dmdmps:>8}'
OUT_PREFETCH = '{pdhps:>8}{pdmps:>8}{pmdhps:>8}{pmdmps:>8}'
OUT_ARC = '{c!s:>10}{p!s:>10}'
ARC_FORMAT = OUT_TOTAL + OUT_DEMAND + OUT_PREFETCH + OUT_ARC

EXTENDEDARC_HEADER = '{anon:^28}{mru:^38}{mfu:^38}{gmru:^38}{gmfu:^38}'
OUT_ANON = '{anons!s:>8}{adevc!s:>10}{amevc!s:>10}'
OUT_MRU = '{mrus!s:>10}{mruh:>8}{mrudevc!s:>10}{mrumevc!s:>10}'
OUT_MFU = '{mfus!s:>10}{mfuh:>8}{mfudevc!s:>10}{mfumevc!s:>10}'
OUT_GMRU = '{gmrus!s:>10}{gmruh:>8}{gmrudevc!s:>10}{gmrumevc!s:>10}'
OUT_GMFU = '{gmfus!s:>10}{gmfuh:>8}{gmfudevc!s:>10}{gmfumevc!s:>10}'
EXTENDEDARC_FORMAT = OUT_ANON + OUT_MRU + OUT_MFU + OUT_GMRU + OUT_GMFU

L2ARC_SIZE = '{l2size!s:>14}{l2hdrsize!s:>14}'
L2ARC_IO = '{l2hps:>8}{l2mps:>8}{l2read!s:>10}{l2write!s:>10}'
L2ARC_FORMAT = L2ARC_SIZE + L2ARC_IO

ZIL_HEADER = '{total:^20}{transactions:^40}{copy:^48}'
ZIL_COMMITS = '{cps:>10}{wps:>10}'
ZIL_ITX = '{itxs:>8}{imnbs!s:>8}{imncs:>8}{imsbs!s:>8}{imscs:>8}'
ZIL_COPIES = '{iibs!s:>8}{iics:>8}{incbs!s:>8}{inccs:>8}{icbs!s:>8}{iccs:>8}'
ZIL_FORMAT = ZIL_COMMITS + ZIL_ITX + ZIL_COPIES

PREFETCH_TOTAL = '{pfhps:>8}{pfmps:>8}'
PREFETCH_INTERN = '{pfchps:>8}{pfcmps:>8}{pfshps:>8}{pfsmps:>8}'
PREFETCH_MISC = '{pfrs:>8}{pfrf:>8}{pfsr:>10}{pfsnr:>12}{pfbs:>10}'
PREFETCH_FORMAT = PREFETCH_TOTAL + PREFETCH_INTERN + PREFETCH_MISC

IO_HEADER = '{bandwidth:^30}{operations:^20}{queue:^32}'
IO_OPS = '{nreads!s:>15}{nwrittens!s:>15}{readss:>10}{writess:>10}'
IO_WAIT = '{awaitq:>8}{waitpct:>8}'
IO_RUN = '{arunq:>8}{runpct:>8}'
IO_FORMAT = IO_OPS + IO_WAIT + IO_RUN

TX_ASSIGN_FORMAT = '{atime:>12}'

def zfsversion():
    cmd = sp.Popen(ZFSVERSION, stdout=sp.PIPE)
    # zfs-0.6.2.19-1.el6.x86_64
    rpm = cmd.communicate()[0]
    ver = rpm.split('-')[1]
    ver = [int(v) for v in ver.split('.')]
    # [0, 6, 2, 19]
    return ver

# Minimal version supportd pool's kstat is zfs-0.6.2.19
def pool_kstat_supported():
    ver = zfsversion()
    if ver[2] < 3:
        try:
            if ver[3] < 19:
                return False
        except IndexError:
            return False
    return True

TXGS_GEN = '{txg:>8}{birth:>20}{state:>8}' if pool_kstat_supported() else '{txg:>8}{state:>8}{birth:>20}'
TXGS_RESRV = '{nreserved:>10}' if pool_kstat_supported() else ''
TXGS_OPS = TXGS_RESRV + '{nread:>10}{nwritten:>12}{reads:>10}{writes:>10}'
TXGS_TIME = '{otime:>12}{qtime:>12}{stime:>12}'
TXGS_FORMAT = TXGS_GEN + TXGS_OPS + TXGS_TIME

IZARCDIR = '/var/log/izarc'

def create_pid():
    if not os.path.exists(IZARCDIR):
        os.mkdir(IZARCDIR)

    PIDFILE = os.path.join(IZARCDIR, str(os.getpid()))
    with open(PIDFILE, 'w'):
        pass

def destroy_pid():
    PIDFILE = os.path.join(IZARCDIR, str(os.getpid()))
    if os.path.exists(PIDFILE):
        os.unlink(PIDFILE)

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

def set_zfsparams(params):
    for item in params:
        with open(os.path.join(MODPARAMPREF, item), 'w') as param:
            param.write(params[item])

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
    KSTAT_NUM_TYPES = 5

    # kstat data types
    KSTAT_DATA_CHAR = 0
    KSTAT_DATA_UINT64 = 4

    def __init__(self, module=None, name=None):
        self._kstat = dict()
        self._raw_data = ''
        if module and name:
            self._file = os.path.join(self.KSTATBASE, module, name)
            self._init_kstat()
        else:
            self._file = ''

    def _init_kstat(self):
        lines = open(self._file).readlines()
        try:
            # Parse common kstat data
            (self._kid, self._type, self._flags, self._ndata, self._data_size,
                self._crtime, self._snaptime) = [int(i, 0) for i in lines[0].split()]
        except IndexError:
            time.sleep(1)
            lines = open(self._file).readlines()
            # Parse common kstat data
            (self._kid, self._type, self._flags, self._ndata, self._data_size,
                self._crtime, self._snaptime) = [int(i, 0) for i in lines[0].split()]

        # Parse the rest of kstat data
        if self._type == self.KSTAT_TYPE_RAW:
            self._init_raw(lines)
        elif self._type == self.KSTAT_TYPE_NAMED:
            self._init_named(lines)
        elif self._type == self.KSTAT_TYPE_INTR:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_TYPE_IO:
            self._init_io(lines)
        elif self._type == self.KSTAT_TYPE_TIMER:
            raise NotImplementedError(self._type)
        elif self._type == self.KSTAT_NUM_TYPES:
            # on zfs < 0.6.2.19 txgs has this type, but now it changed to KSTAT_TYPE_RAW
            self._init_raw(lines)
        else:
            raise NotImplementedError(self._type)

    def _init_named(self, lines):
        for line in lines[2:]:
            try:
                name, data, value = line.split()
            except ValueError:
                name, metric, data, value = line.split()
                name += ' ' + metric
            data = int(data)
            if data == self.KSTAT_DATA_CHAR:
                pass
            elif data == self.KSTAT_DATA_UINT64:
                self._kstat[name] = int(value, 0)
            else:
                pass

    def _init_io(self, lines):
        headers = lines[1].split()
        data = [long(d) for d in lines[2].split()]
        self._kstat = dict(zip(headers, data))

    def _init_raw(self, lines):
        self._raw_data = lines[1:]
        return self._raw_data

    def __sub__(self, other):
        if not isinstance(other, kstat):
            return NotImplemented
        diff = deepcopy(self)
        diff._snaptime = self._snaptime - other._snaptime
        for item in self._kstat:
            diff._kstat[item] = Integer(self._kstat[item] - other._kstat[item])
        return diff

    def __str__(self):
        text = 'kid {self._kid}\n'.format(self=self)
        text = pp.pformat(self._kstat)
        return text


class tx_assign(kstat):
    MAX_TIME = 40
    TX_ASSIGN_FORMAT = ''

    def __init__(self, pool):
        self._pool = pool
        if pool_kstat_supported():
            super(tx_assign, self).__init__(os.path.join('zfs', self._pool), 'dmu_tx_assign')
        else:
            super(tx_assign, self).__init__(os.path.join('zfs'), 'dmu_tx_assign-{0}'.format(self._pool))

    @classmethod
    def acronym(cls):
        out = '\n----- TX_ASSIGN Acronym -----\n'
        out += 'atime/s     - Average time (in ns) that transactions waited until enter into the new txg'
        return out

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._kstat.copy()
        all_tx = 0
        tx_time = 0
        for i in range(self.MAX_TIME):
            units = 'ns' if pool_kstat_supported() else 'us'
            spent_time = pow(2, i)
            t = '{0} {1}'.format(spent_time, units)
            raw[t] = self._kstat.get(t, 0) / delta
            all_tx += raw[t]
            tx_time += (spent_time * raw[t])

        try:
            raw['atime'] = float(tx_time) / float(all_tx)
        except ZeroDivisionError:
            raw['atime'] = 0
        return raw

    def headers(self):
        header = '\n'
        header += TX_ASSIGN_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return TX_ASSIGN_FORMAT.format(**raw)


class txgs(kstat):
    def __init__(self, pool):
        self._pool = pool
        if pool_kstat_supported():
            super(txgs, self).__init__(os.path.join('zfs', self._pool), 'txgs')
        else:
            super(txgs, self).__init__(os.path.join('zfs'), 'txgs-{0}'.format(self._pool))

    @classmethod
    def acronym(cls):
        out = '\n----- TXGS Acronym -----\n'
        out += 'txg           - Transaction group id\n'
        out += 'birth         - Transaction group birth timestamp in ns\n'
        out += 'state         - Transaction group state (Open/Quiescing/Syncing)\n'
        out += 'nread         - Number of bytes read\n'
        out += 'nwritten      - Number of bytes written\n'
        out += 'reads         - Number of read operations\n'
        out += 'writes        - Number of write operations\n'
        out += 'otime         - Time spent in Open state in ns\n'
        out += 'qtime         - Time spent in Quiescing state in ns\n'
        out += 'wtime         - Time spent in Wait (between Quiescing and Syncing) state in ns\n'
        out += 'stime         - Time spent in Syncing state in ns\n'
        return out

    def data(self):
        return self._raw_data

    def latest_data(self):
        data = self._raw_data[0] if len(self._raw_data) > 5 else ''
        for line in self._raw_data[-5:]:
            data += line
        return data

    def __str__(self):
        data = ''
        for line in self._raw_data:
            data += line
        return data


class io(kstat):
    def __init__(self, pool):
        self._pool = pool
        super(io, self).__init__(os.path.join('zfs', self._pool), 'io')

    @classmethod
    def acronym(cls):
        out = '\n----- IO Acronym -----\n'
        out += 'read/s          - Number of bytes/operations read per second\n'
        out += 'write/s         - Number of bytes/operations written per second\n'
        out += 'awaitq          - Average length of wait queue\n'
        out += 'waitq%          - Percentage of time spent in wait queue at given time interval\n'
        out += 'arunq           - Average length of run queue\n'
        out += 'runq%           - Percentage of time spent in run queue at given time interval\n'
        # out += 'wtime           - cumulative wait (pre-service) time (sec)\n'
        # out += 'wlentime        - cumulative wait length*time product (sec)\n'
        # out += 'wlastupdate     - last time wait queue changed (sec)\n'
        # out += 'rtime           - cumulative run (service) time (sec)\n'
        # out += 'rlentime        - cumulative run length*time produc (sec)\n'
        # out += 'rlastupdate     - last time run queue changed (sec)\n'
        # out += 'wcnt            - count of elements in wait state\n'
        # out += 'rcnt            - count of elements in run state\n'
        return out

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._kstat.copy()
        raw['nreads'] = Integer(self._kstat['nread'] / delta)
        raw['nwrittens'] = Integer(self._kstat['nwritten'] / delta)
        raw['readss'] = self._kstat['reads'] / delta
        raw['writess'] = self._kstat['writes'] / delta
        # Accumulated time and queue length statistics.
        #
        # Accumulated time statistics are kept as a running sum
        # of "active" time.  Queue length statistics are kept as a
        # running sum of the product of queue length and elapsed time
        # at that length
        # At each change of state (entry or exit from the queue),
        # we add the elapsed time (since the previous state change)
        # to the active time if the queue length was non-zero during
        # that interval; and we add the product of the elapsed time
        # times the queue length to the running length*time sum.
        #
        # All times are 64-bit nanoseconds (hrtime_t), as returned by gethrtime().
        # raw['wtime'] = long(self._kstat['wtime'])
        # raw['wlentime'] = int(self._kstat['wlentime'])
        # raw['wupdate'] = int(self._kstat['wupdate'])
        # raw['rtime'] = int(self._kstat['rtime'])
        # raw['rlentime'] = int(self._kstat['rlentime'])
        # raw['rupdate'] = int(self._kstat['rupdate'])
        # raw['wcnt'] = self._kstat['wcnt']
        # raw['rcnt'] = self._kstat['rcnt']

        getcontext().prec = 3
        try:
            raw['awaitq'] = Decimal(self._kstat['wlentime']) / Decimal(self._kstat['wtime'])
        except (InvalidOperation, DivisionByZero):
            raw['awaitq'] = 0
        raw['waitpct'] = Decimal(self._kstat['wtime']) * 100 / Decimal(self._snaptime)
        try:
            raw['arunq'] = Decimal(self._kstat['rlentime']) / Decimal(self._kstat['rtime'])
        except (InvalidOperation, DivisionByZero):
            raw['arunq'] = 0
        raw['runpct'] = Decimal(self._kstat['rtime']) * 100 / Decimal(self._snaptime)
        return raw

    def headers(self):
        header = '\n'
        header += IO_HEADER.format(**HEADER_NAMES)
        header += '\n'
        header += IO_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return IO_FORMAT.format(**raw)


class prefetch(kstat):
    def __init__(self):
        super(prefetch, self).__init__('zfs', 'zfetchstats')

    @classmethod
    def acronym(cls):
        out = '\n----- Prefetch Acronym -----\n'
        out += 'hit/s      - Hits per second\n'
        out += 'miss/s     - Misses per second\n'
        out += 'chit/s     - Colinear hits per second\n'
        out += 'cmiss/s    - Colinear misses per second\n'
        out += 'shit/s     - Stride hits per second\n'
        out += 'smiss/s    - Stride misses per second\n'
        out += 'rsuc/s     - Reclaim successes per second\n'
        out += 'rfail/s    - Reclaim failures per second\n'
        out += 'sbogus/s   - Bogus streams per second\n'
        out += 'sreset/s   - Streams resets per second\n'
        out += 'snoreset/s - Streams noresets per second\n'
        return out

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._kstat.copy()
        raw['pfhps'] = self._kstat['hits'] / delta
        raw['pfmps'] = self._kstat['misses'] / delta
        raw['pfchps'] = self._kstat['colinear_hits'] / delta
        raw['pfcmps'] = self._kstat['colinear_misses'] / delta
        raw['pfshps'] = self._kstat['stride_hits'] / delta
        raw['pfsmps'] = self._kstat['stride_misses'] / delta
        raw['pfrs'] = self._kstat['reclaim_successes'] / delta
        raw['pfrf'] = self._kstat['reclaim_failures'] / delta
        raw['pfsr'] = self._kstat['streams_resets'] / delta
        raw['pfsnr'] = self._kstat['streams_noresets'] / delta
        raw['pfbs'] = self._kstat['bogus_streams'] / delta
        return raw

    def headers(self):
        header = '\n'
        header += PREFETCH_FORMAT.format(**METRIC_NAMES)
        return header

    def data(self):
        raw = self.compute()
        return PREFETCH_FORMAT.format(**raw)


class arcstats(kstat):
    def __init__(self):
        super(arcstats, self).__init__('zfs', 'arcstats')
        self._arcstats = dict()
        for name in self._kstat:
            self._arcstats[name] = Integer(self._kstat[name])

    @classmethod
    def acronym(cls):
        out = '\n----- Arcstats Acronym -----\n'
        out += 'hit/s      - Hits per second\n'
        out += 'miss/s     - Misses per second\n'
        out += 'datah/s    - Data hits per second\n'
        out += 'datam/s    - Data misses per second\n'
        out += 'metah/s    - Metadata hits per second\n'
        out += 'metam/s    - Metadata misses per second\n'
        out += 'deltac     - Delta of target ARC size at given interval\n'
        out += 'deltap     - Delta of target MRU size at given interval\n'
        return out

    def __sub__(self, other):
        if not isinstance(other, arcstats):
            return NotImplemented
        diff = deepcopy(self)
        diff._snaptime = self._snaptime - other._snaptime
        for item in self._arcstats:
            diff._arcstats[item] = Integer(self._arcstats[item] - other._arcstats[item])
        return diff

    def __str__(self):
        return str(self._arcstats)

    def __repr__(self):
        return str(self._arcstats)

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
    @classmethod
    def acronym(cls):
        out = '\n----- Extended Arcstats Acronym -----\n'
        out += 'hit/s      - Hits per second\n'
        out += 'mevict     - Evictable metadata\n'
        out += 'devict     - Evictable data\n'
        out += 'dltsize    - Delta of size at given interval\n'
        return out

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._arcstats.copy()
        raw['anons'] = self._arcstats['anon_size']
        raw['adevc'] = Integer(self._arcstats['anon_evict_data'])
        raw['amevc'] = Integer(self._arcstats['anon_evict_metadata'])
        raw['mrus'] = self._arcstats['mru_size']
        raw['mruh'] = self._arcstats['mru_hits'] / delta
        raw['mrudevc'] = Integer(self._arcstats['mru_evict_data'])
        raw['mrumevc'] = Integer(self._arcstats['mru_evict_metadata'])
        raw['mfus'] = self._arcstats['mfu_size']
        raw['mfuh'] = self._arcstats['mfu_hits'] / delta
        raw['mfudevc'] = Integer(self._arcstats['mfu_evict_data'])
        raw['mfumevc'] = Integer(self._arcstats['mfu_evict_metadata'])
        raw['gmrus'] = self._arcstats['mru_ghost_size']
        raw['gmruh'] = self._arcstats['mru_ghost_hits'] / delta
        raw['gmrudevc'] = Integer(self._arcstats['mru_ghost_evict_data'])
        raw['gmrumevc'] = Integer(self._arcstats['mru_ghost_evict_metadata'])
        raw['gmfus'] = self._arcstats['mfu_ghost_size']
        raw['gmfuh'] = self._arcstats['mfu_ghost_hits'] / delta
        raw['gmfudevc'] = Integer(self._arcstats['mfu_ghost_evict_data'])
        raw['gmfumevc'] = Integer(self._arcstats['mfu_ghost_evict_metadata'])
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
    @classmethod
    def acronym(cls):
        out = '\n----- L2ARC Acronym -----\n'
        out += 'hit/s          - Hits per second\n'
        out += 'miss/s         - Misses per second\n'
        out += 'read/s         - Read bytes per second\n'
        out += 'write/s        - Write bytes per second\n'
        out += 'dltsize        - Delta of L2ARC size at given interval\n'
        out += 'dlthdrsize     - Delta of L2ARC header size (ARC overhead) at given interval\n'
        return out

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._arcstats.copy()
        raw['l2size'] = self._arcstats['l2_size']
        # L2ARC headers (ARC overhead)
        raw['l2hdrsize'] = self._arcstats['l2_hdr_size']
        raw['l2hps'] = self._arcstats['l2_hits'] / delta
        raw['l2mps'] = self._arcstats['l2_misses'] / delta
        raw['l2read'] = Integer(self._arcstats['l2_read_bytes'] / delta)
        raw['l2write'] = Integer(self._arcstats['l2_write_bytes'] / delta)
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

    @classmethod
    def acronym(cls):
        out = '\n----- ZIL Acronym -----\n'
        out += 'commit/s   - Commit requests per second\n'
        out += 'wrcomm/s   - Number of times the ZIL has been flushed to stable storage per second\n'
        out += 'itx/s      - Number of transactions (reads, writes, renames, etc.) that have been commited per second\n'
        out += '*** Transactions which have been allocated to the "normal" (i.e. not slog) storage pool ***\n'
        out += 'imsnb/s    - Accumulates the actual log record sizes per second\n'
        out += 'imsnc/s    - Writes counter per second\n'
        out += '*** Transactions which have been allocated to the "slog" storage pool ***\n'
        out += 'imssb/s    - Accumulates the actual log record sizes per second per second\n'
        out += 'imssc/s    - Writes counter per second\n'
        out += '*** Counters for three different ways of writes handling ***\n'
        out += 'iib/s      - Accumulates the length of the indirect data, not the actual log record sizes per second\n'
        out += 'iic/s      - Indirect data writes counter per second\n'
        out += 'incb/s     - Accumulates the length of the needcopy data, not the actual log record sizes per second\n'
        out += 'incc/s     - Needcopy data writes counter per second\n'
        out += 'icb/s      - Accumulates the length of the copied data, not the actual log record sizes per second\n'
        out += 'icc/s      - Copied data writes counter per second\n'
        return out

    def compute(self):
        delta = self._snaptime / NANOSEC
        raw = self._kstat.copy()
        # Number of times a ZIL commit (e.g. fsync) has been requested.
        raw['cps'] = self._kstat['zil_commit_count'] / delta
        # Number of times the ZIL has been flushed to stable storage.
        # This is less than zil_commit_count when commits are "merged".
        raw['wps'] = self._kstat['zil_commit_writer_count'] / delta
        # Number of transactions (reads, writes, renames, etc.) that have been commited.
        raw['itxs'] = self._kstat['zil_itx_count'] / delta
        # Writes are handled in three different ways:
        #
        # WR_INDIRECT:
        #    In this mode, if we need to commit the write later, then the block
        #    is immediately written into the file system (using dmu_sync),
        #    and a pointer to the block is put into the log record.
        #    When the txg commits the block is linked in.
        #    This saves additionally writing the data into the log record.
        #    There are a few requirements for this to occur:
        #  - write is greater than zfs/zvol_immediate_write_sz
        #  - not using slogs (as slogs are assumed to always be faster
        #    than writing into the main pool)
        #  - the write occupies only one block
        # WR_COPIED:
        #    If we know we'll immediately be committing the
        #    transaction (FSYNC or FDSYNC), the we allocate a larger
        #    log record here for the data and copy the data in.
        # WR_NEED_COPY:
        #    Otherwise we don't allocate a buffer, and *if* we need to
        #    flush the write later then a buffer is allocated and
        #    we retrieve the data using the dmu.

        # Note that "bytes" accumulates the length of the transactions
        # (i.e. data), not the actual log record sizes.
        raw['iibs'] = Integer(self._kstat['zil_itx_indirect_bytes'] / delta)
        raw['iics'] = self._kstat['zil_itx_indirect_count'] / delta
        raw['incbs'] = Integer(self._kstat['zil_itx_needcopy_bytes'] / delta)
        raw['inccs'] = self._kstat['zil_itx_needcopy_count'] / delta
        raw['icbs'] = Integer(self._kstat['zil_itx_copied_bytes'] / delta)
        raw['iccs'] = self._kstat['zil_itx_copied_count'] / delta
        # Transactions which have been allocated to the "normal"
        # (i.e. not slog) storage pool. Note that "bytes" accumulate
        # the actual log record sizes - which do not include the actual
        # data in case of indirect writes.
        raw['imnbs'] = Integer(self._kstat['zil_itx_metaslab_normal_bytes'] / delta)
        raw['imncs'] = self._kstat['zil_itx_metaslab_normal_count'] / delta
        # Transactions which have been allocated to the "slog" storage pool.
        # If there are no separate log devices, this is the same as the
        # "normal" pool.
        raw['imsbs'] = Integer(self._kstat['zil_itx_metaslab_slog_bytes'] / delta)
        raw['imscs'] = self._kstat['zil_itx_metaslab_slog_count'] / delta
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

def cycle(stats, interval, count, verbose, debug, timestamp, **kwargs):
    if verbose:
        print stats.acronym()
    step = 1 if count else 0
    count = count if count else 1
    cur = stats(**kwargs)
    if debug:
        print cur
    time.sleep(1)
    lines = 0
    while count:
        count -= step
        prev = cur
        cur = stats(**kwargs)
        delta = cur - prev
        if lines % 20 == 0:
            header = delta.headers()
            if header:
                print header

        probetime = '    {0}'.format(datetime.now()) if timestamp else ''
        data = delta.data()
        if data:
            lines += 1
            print data + probetime
        if count:
            time.sleep(interval)

def raw_history(stats, interval, verbose, timestamp, **kwargs):
    if verbose:
        print stats.acronym()
    time.sleep(5)
    cur = stats(**kwargs)
    print cur
    while True:
        cur = stats(**kwargs)
        probetime = '    {0}\n'.format(datetime.now()) if timestamp else ''
        print probetime + cur.latest_data()
        time.sleep(interval)

def execute(args):
    if args.outfile:
        try:
            out = open(args.outfile, "w")
            sys.stdout = out
        except IOError:
            sys.stderr.write("Cannot open %s for writing\n" % args.outfile)
            sys.exit(1)

    if args.summary:
        if args.time:
            print datetime.now()
        print arcstats().summary()
    elif args.parm:
        print zfsparams()
    elif args.zil:
        cycle(zil, args.interval, args.count, args.verbose, args.debug, args.time)
    elif args.l2arc:
        cycle(l2arc, args.interval, args.count, args.verbose, args.debug, args.time)
    elif args.extend:
        cycle(extendarc, args.interval, args.count, args.verbose, args.debug, args.time)
    elif args.prefetch:
        cycle(prefetch, args.interval, args.count, args.verbose, args.debug, args.time)
    elif args.pool:
        if pool_kstat_supported():
            set_zfsparams(dict(zfs_read_history='100',
                               zfs_read_history_hits='1',
                               zfs_txg_history='100'))
        if args.io:
            if pool_kstat_supported():
                cycle(io, args.interval, args.count, args.verbose, args.debug, args.time, pool=args.pool)
            else:
                print "--io option supported for zfs >= 0.6.2.19 only"
        elif args.tx_assign:
            cycle(tx_assign, args.interval, args.count, args.verbose, args.debug, args.time, pool=args.pool)
        elif args.txgs:
            raw_history(txgs, args.interval, args.verbose, args.time, pool=args.pool)
        elif args.read:
            if pool_kstat_supported():
                pass
            else:
                print "--read option supported for zfs >= 0.6.2.19 only"
    # elif args.debug:
    #     print arcstats()
    #     as1 = arcstats()
    #     time.sleep(1)
    #     as2 = arcstats()
    #     print (as2 - as1)
    else:
        cycle(arcstats, args.interval, args.count, args.verbose, args.debug, args.time)


def main():
    parser = argparse.ArgumentParser(
        description='Report various ZFS ARC statistics')
    parser.add_argument('-v', '--version', action='version', version=VERSION)
    parser.add_argument('-d', '--debug',
        help='help debug this tool', action='store_true')
    parser.add_argument('--verbose',
        help='Statistics acronym', action='store_true')
    parser.add_argument('-o', '--outfile',
        help='Redirect output to the specified file', nargs='?')
    parser.add_argument('-t', '--time',
        help='print timestamp', action='store_true')
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
    parser.add_argument('-f', '--prefetch',
        help='print prefetch statistics', action='store_true')
    parser.add_argument('--pool',
        help='pool name to report statistics', nargs='?')
    parser.add_argument('--tx_assign',
        help='print POOL\'s dmu_tx_assign statistics', action='store_true')
    parser.add_argument('--txgs',
        help='print POOL\'s txgs statistics', action='store_true')
    parser.add_argument('--io',
        help='print POOL\'s io statistics', action='store_true')
    # parser.add_argument('--read',
    #     help='print POOL\'s read statistics', action='store_true')

    parser.add_argument('interval',
        help='seconds between probes', type=int, nargs='?', default=1)
    parser.add_argument('count',
        help='number of probes to run', type=int, nargs='?', default=0)

    args = parser.parse_args()
    create_pid()
    execute(args)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass

    destroy_pid()
    pids = glob(os.path.join(IZARCDIR, '*'))
    if pool_kstat_supported() and not pids:
        set_zfsparams(dict(zfs_read_history='0',
                           zfs_read_history_hits='0',
                           zfs_txg_history='0'))
