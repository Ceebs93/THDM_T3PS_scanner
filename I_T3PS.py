#!/usr/bin/env python2
# coding=utf-8
#
# Copyright 2014-2016 Vinzenz Maurer
# Edited by David Englert and Ciara Byers
#
# This file is part of T3PS.
#
# T3PS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# T3PS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with T3PS.  If not, see <http://www.gnu.org/licenses/>.
"""
Parameter scanner using the multiprocessing module.

Run this file the command line argument `--help' for details on its usage or
consult the manual.
"""

# Ensuring division of two integers returns a float, i.e. 1/2 -> 0.5
from __future__ import division

__program_name__ = "T3PS"
__major_version__ = 1
__version__ = "1.0.2"

import os
import sys

if sys.platform == 'win32':
    # win32 does not support os.fork(), making multiprocessing much more
    #   difficult to use
    # in particular, point processors would have to be initialized multiple
    #   times or have to save their own derived config in picklable form
    sys.exit("Operating system not supported")

# turning off output buffering
# (important for progress bars and general responsiveness)
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

# ---------------------------------------------------------------------- #
# | For easier distribution, this program is one monolithic file.      | #
# | Thus, to compensate for missing structure we include a searchable  | #
# ------------------------ Table of Contents: -------------------------- #
# | [A] Imports
# | [B] Constants
# | [C] Command Line Interface
# | [D] Classes, Methods and Functions
# |  [D.01] Contexts
# |  [D.02] Utility
# |  [D.03] Progress Display
# |  [D.04] Iterator Handling
# |  [D.05] File Input/Output
# |  [D.06] Configuration Input (Parameter Range, Formulas, Preprocessing)
# |  [D.07] Randomness
# |  [D.08] Markov Chain Monte Carlo Functions
# |  [D.09] Explorer Mode Functions
# |  [D.10] Point Processing
# |  [D.11] Program Flow Control
# | [E] Configuration File Handling
# |  [E.1] Processing Directive
# |  [E.1.1] General Directives
# |  [E.1.2] Shared Directives
# |  [E.1.3] Specific Directives
# |  [E.2] Review by user
# | [F] Manager Mode Setup
# | [G] Scanning Code
# |  [G.1] Test Mode
# |  [G.2] Scan Mode
# |  [G.3] MCMC Mode
# |  [G.4] Optimize Mode
# |  [G.5] Explorer Mode
# --------------------------------------------------------------------- #
# | All parts can be searched and found by their [label].             | #
# --------------------------------------------------------------------- #

# --------------- #
# | [A] Imports | #
# --------------- #

import array
import argparse
import atexit
import ConfigParser
import cPickle as pickle

# Parameter ranges will use Decimals instead of floats to minimize floating
# point misbehaviour in sums
#   e.g. 1.3 + 0.1 == 1.4000000000000001, 0.7 + 0.1 == 0.7999999999999999
import decimal
import gc
import imp
import hashlib
import itertools
import keyword
import math
import mmap
import multiprocessing
import os.path
import random
import re
import shutil
import signal
import string
import StringIO  # cannot use cStringIO, as we change attributes in it later
import textwrap
import time
import traceback
import warnings

# for determining terminal size
import fcntl
import struct
import termios

from collections import namedtuple, Sequence, deque
from contextlib import contextmanager
from decimal import Decimal
from glob import glob
from multiprocessing.managers import BaseManager
from tempfile import mkdtemp


# ----------------- #
# | [B] Constants | #
# ----------------- #

MODE_SCAN = 0
MODE_MCMC = 1
MODE_OPTIMIZE = 2
MODE_EXPLORER = 3
MODE_TEST = 1001

PARSPACE_MODE_GRID = 10
PARSPACE_MODE_FILE = 11
PARSPACE_MODE_SCATTER = 12

DEFAULT_PORT = 31415

mode_names = {
    "scan": MODE_SCAN,
    "mcmc": MODE_MCMC,
    "test": MODE_TEST,
    "explorer": MODE_EXPLORER,
    "optimize": MODE_OPTIMIZE
}
parspace_mode_names = {
    "grid": PARSPACE_MODE_GRID,
    "file": PARSPACE_MODE_FILE,
    "scatter": PARSPACE_MODE_SCATTER,
}


# If the first column contains this character it is an invalid line
ERROR_MARKER = "E"

# ------------------------------ #
# | [C] Command Line Interface | #
# ------------------------------ #

# always listen to SIGINT!
signal.signal(signal.SIGINT, signal.default_int_handler)

parser = argparse.ArgumentParser(
    description='Perform parameter scans over multiple processors ' +
    '(local and remote).'
)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%s v%s' % (__program_name__, __version__)
)
parser.add_argument(
    'input_file', type=argparse.FileType("r"), nargs='*',
    help='configuration file specifying the parameter scan'
)
mode_choices = [
    k for k, v in sorted(
        mode_names.iteritems(),
        key=lambda k_and_v: k_and_v[1]
    )
]
parser.add_argument(
    '--mode',
    type=str, action='store',
    choices=mode_choices,
    help='overwrite the mode specified in the configuration file'
)
parser.add_argument(
    '-o', '--output_dir', action='store',
    help='write all output to this directory (will be created if needed)'
)
parser.add_argument(
    '-p', '--port', type=int, action='store',
    help='listen on this port (worker mode)'
)
parser.add_argument(
    '--pars', type=float, action='append', nargs='+', metavar="VALUE",
    help='process the points given by these parameter values, ' +
    'can be given multiple times (test mode)'
)
parser.add_argument(
    '-P', '--profiling', action='store_true',
    help='enable profiling info output (test mode)'
)
parser.add_argument(
    '--experimental', action="store_true", help=argparse.SUPPRESS
)
parser.add_argument(
    '-D', '--debug', action="store_true", help=argparse.SUPPRESS
)
parser.add_argument(
    '--randomseed', type=int, action='store',
    help=argparse.SUPPRESS
)

# Considerations on determinism and randomness
# The aim of this tool is that given a specific random seed, the exact same
# behavior should be exhibited.
# Thus the following "arbitrary" behavior must be excised:
#   - set iteration order
#   - ...

cli_arguments = parser.parse_args()

print "#", __program_name__, "v%s" % __version__

# -------------------------------------- #
# | [D] Classes, Methods and Functions | #
# -------------------------------------- #

# -------------------- #
# |  [D.01] Contexts | #
# -------------------- #

# Context manager facilitates code such as:
#   <code executing in directory 1>
#   with pushd(dir2):
#       <code executing in directory dir2>
#   <code executing in directory 1>
# advantages: exception safe handling


@contextmanager
def pushd(newdir):
    """
    Context that causes temporary change to another directory.
    
    Parameters
    ----------
    newdir : string
        Path to a directory we wish to change to temporarily

    Returns
    -------
    None
    """

    currentdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(currentdir)


@contextmanager
def TemporaryDirectory(keep_dir=False):
    """
    Context that creates a temporary directory during live block execution.

    Parameters
    ----------
    keep_dir: Boolean 
        Indicates whether or not to delete temp directory immediately after
        leaving block, occurs only at the end of function's execution

    Returns
    -------
    None
    """

    tmpdir = mkdtemp()
    if keep_dir:
        atexit.register(shutil.rmtree, tmpdir)
    try:
        yield tmpdir
    finally:
        if not keep_dir:
            shutil.rmtree(tmpdir)


@contextmanager
def ConfirmedExitOnInterrupt():
    """
    Context that makes a block require two SIGINTs within 5 seconds to stop.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    received_signals = []
    # Setting default functionality of signal.get
    old_handler = signal.getsignal(signal.SIGINT)
    if old_handler == signal.SIG_DFL:
        old_handler = signal.default_int_handler

    def handler(sig, frame, received_signals=received_signals,
                old_handler=old_handler):
        # only abort on two interrupts
        if received_signals and time.time() - max(received_signals) < 5:
            old_handler(sig, frame)
        received_signals.append(time.time())

    # replace the handler with ours
    signal.signal(signal.SIGINT, handler)

    try:
        yield  # execute the with block
    except KeyboardInterrupt:
        print  # put those ^C on their own lines
        exit_program()
    finally:
        # restore the old handler
        signal.signal(signal.SIGINT, old_handler)


@contextmanager
def DelayedKeyboardInterrupt(wait_message=""):
    """Context making a block atomic with respect to interrupts.

    It does this by delaying it until after block execution.
    """
    # NOTE: it is not necessary to also intercept EOFError or similar, since
    #   those only come with raw_input

    # make a list, such that referencing it in the handler works
    received_signal = [False]
    old_handler = signal.getsignal(signal.SIGINT)

    def handler(sig, frame, received_signal=received_signal,
                wait_message=wait_message):
        # save the signal and instead show a message
        received_signal[0] = (sig, frame)
        if wait_message:
            print wait_message

    # replace the handler with ours
    signal.signal(signal.SIGINT, handler)

    try:
        yield  # execute the with block
    finally:
        # restore the old handler and let it its job if needed
        signal.signal(signal.SIGINT, old_handler)
    if received_signal[0]:
        old_handler(*received_signal[0])


class TimeoutError(Exception):

    """Custom exception class signaling a run out timeout."""

    pass


@contextmanager
def TimeLimit(seconds):
    """Context for blocks that are aborted if they take too long."""
    def handler(*_args):
        raise TimeoutError("Time limit exceeded")

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(int(round(seconds)))
    try:
        yield
    finally:
        signal.alarm(0)


# ------------------- #
# |  [D.02] Utility | #
# ------------------- #

def prod(iterable):
    """Calculates the Product of a number, i.e. the analogue of the sum
    function but for multiplication.
    
    Parameters 
    ----------
    iterable : some iterable object such as a list or array containing numbers

    Returns
    -------
    r : float/integer
        returns the product of the input which will have the same typing as the
        values in the input iterable object.
    """
    r = 1
    for x in iterable:
        r *= x
    return r


def list_sum(lists):
    """
    Concatenate lists together, similar to sum function.
    
    Parameters
    ----------
    lists : a list/array of lists
        each element should be a list that the user wishes to concatenate onto
        the end of the previous element.
"""
    # For reference (cannot use named argument):
    #   sum(items, start_value)
    return sum(lists, [])


def inverseerf(x, sqrt=math.sqrt, pi=math.pi, log=math.log, a=0.147):
    """
    Inverse erf function.

    Approximation taken from Winitzki, Sergei (6 February 2008)
    "A handy approximation for the error function and its inverse" (PDF).
    http://sites.google.com/site/winitzki/sergei-winitzkis-files/erf-approx.pdf
    This has relative precision better than 2e-3 for x in (-1, 1).

    Parameters
    ----------
    x : float

    sqrt : func
        default = math.sqrt, can provide alternate function for sqrt if desired

    pi : float
        default = math.pi, can provide alternative value for pi if desired

    log : func
        default = math.log, can provide alternative function for log if desired

    a : float
        default = 0.147, can provide alternative value if desired

    Returns
    -------
    math.copysign(..) : float
        value of sqrt action performed with the same sign as x assigned to it
    """
    return math.copysign(
        sqrt(
            sqrt((2 / (pi * a) + log(1 - x * x) / 2) ** 2 - log(1 - x * x) / a)
            - (2 / (pi * a) + log(1 - x * x) / 2)
        ),
        x
    )


# thus this has accuracy 0.3%
def inversenormalcdf(p, sqrt=math.sqrt):
    """
    Inverse cumulative distribution function for std. normal distribution.

    Has about the same level of accuracy as the InverseCDF function in
    Wolfram Mathematica
    """
    return sqrt(2) * inverseerf(2 * p - 1)


def maybefloat(s):
    """Convert to float if possible, otherwise to string."""
    try:
        return float(s)
    except Exception:
        return str(s)


def maybedecimal(s):
    """Convert to Decimal if possible, otherwise to string."""
    try:
        return Decimal(s)
    except Exception:
        return str(s)


def is_number(s):
    """Check if string is valid decimal number."""
    try:
        Decimal(s)
        return True
    except decimal.InvalidOperation:
        return False
    except ValueError:
        return False


def print_table(names, values):
    """Print named values in nice table form to console."""
    def to_str(v):
        return repr(v) if isinstance(v, float) else str(v)

    max_data_length = max(max(len(to_str(x)) for x in values), 3)
    max_name_length = max(max(len(x) for x in names), 3)
    max_width = terminal_size()[0]

    # +3 for " | " and +2 for ": "
    headered_value_length = max_data_length + max_name_length + 3 + 2
    # a // b: integer division (rounds down)
    row_length = max(1, max_width // headered_value_length)

    padded_length = max(len(names), len(values))
    row_length = min(row_length, padded_length)
    if padded_length % row_length != 0:
        # print full rows
        padded_length += row_length - (padded_length % row_length)

    values = [to_str(v) for v in values]
    padded_data = list(pad_list(values, "---", padded_length))
    padded_names = list(pad_list(names, "---", padded_length))

    # There is
    #   ": "(2 chars) for each column -> weight row_length
    #   " | "(3 chars) between each column -> weight row_length-1
    hline = "-" * (
        (max_data_length + max_name_length) * row_length +
        2 * row_length +
        3 * (row_length - 1)
    )

    print hline
    # go from 0 to padded_length (exclusive) with step size row_length
    for i in range(0, padded_length, row_length):
        print " | ".join(
            padded_names[i + j].ljust(max_name_length) + ": " +
            padded_data[i + j].rjust(max_data_length)
            for j in range(row_length)
        )
    print hline


def split_by_comma(s):
    """Split by base level comma.

    That means not by comma inside brackets or strings.
    """
    # keep track in what string type we are in currently and what kind of
    #   brackets are open
    string_type = None
    brackets = [0, 0, 0]  # (), [], {}

    token = ""
    parts = []
    for c in s:
        if c == ',':
            # only split by lowest level commas that are not in strings
            if not any(brackets) and string_type is None:
                parts.append(token)
                token = ""
                continue
        elif c == string_type:
            string_type = None
        elif c in ("'", '"'):
            if string_type is None:
                string_type = c
        elif not string_type:
            if c == "(":
                brackets[0] += 1
            elif c == "[":
                brackets[1] += 1
            elif c == "{":
                brackets[2] += 1
            elif c == ")":
                brackets[0] -= 1
            elif c == "]":
                brackets[1] -= 1
            elif c == "}":
                brackets[2] -= 1
        token += c
    parts.append(token)
    return parts


# ---------------------------- #
# |  [D.03] Progress Display | #
# ---------------------------- #

def terminal_size():
    """Determine terminal size as (width, height)."""
    if not sys.stdout.isatty():
        return (80, 24)

    height, width, _height_in_pixels, _width_in_pixels = struct.unpack(
        'HHHH',
        fcntl.ioctl(
            1,  # fd == 1 -> stdout
            termios.TIOCGWINSZ, struct.pack('HHHH', 0, 0, 0, 0)
        )
    )
    return width, height


def horizontal_line():
    """Print horizontal line spanning full width of terminal."""
    return "-" * terminal_size()[0]


def draw_progress(current, maximum, info="", width=None, abortable=True):
    """Draw progress bar to terminal with additional informational message.

    Also includes message if current action is abortable.
    """
    if not sys.stdout.isatty():
        # there's no use to show anything if there's no one looking
        return
    if not hasattr(draw_progress, "last_progress"):
        draw_progress.last_progress = None
        draw_progress.last_progress_time = time.time()
        draw_progress.last_info_length = 0
    # only draw once every 0.1 seconds or when we reached 100%
    elif time.time() - draw_progress.last_progress_time < 0.1 and \
            current < maximum:
        return

    draw_progress.last_progress_time = time.time()

    if width is None:
        width = terminal_size()[0]
    progress = current / (1.0 * maximum)

    # make display of current as long as maximum
    maximum = str(maximum)
    current = str(current).rjust(len(maximum))
    # and leave enough space for "100.00%" (7 chars)
    progress_info = "({}/{}) ~ {:>7.2%}".format(current, maximum, progress)
    # this way the bar never has to be resized

    # -2 for [], -1 for space between bar and info, -10 courtesy
    bar_length = width - 2 - 1 - len(progress_info) - 10

    num_segments = min(int(bar_length * progress), bar_length)
    wrapped_lines = []
    if info:
        for line in info.split("\n"):
            if len(line) <= width:
                wrapped_lines.append(line)
            else:
                wrapped_lines.extend(
                    textwrap.wrap(line, width=width, replace_whitespace=False)
                )
    wrapped_lines.append(
        "current scan: %s" % abbreviate_path(
            draw_progress.scan_file,
            width=width - len("current scan: ")
        )
    )
    if abortable:
        wrapped_lines.append("(press ctrl+c twice within 5 sec to abort)")

    wrapped_lines_length = len(wrapped_lines)
    wrapped_lines = "\n".join(wrapped_lines)

    if draw_progress.last_progress is None or \
       draw_progress.last_progress > progress or \
       draw_progress.last_max != maximum:
        # setup full progress bar space if there is no sufficient previous
        #   progress bar shown at the moment
        for _ in range(1 + wrapped_lines_length):
            print ""
    elif draw_progress.last_info_length < wrapped_lines_length:
        # only setup missing space
        for _ in range(wrapped_lines_length - draw_progress.last_info_length):
            print ""

    draw_progress.last_progress = progress
    draw_progress.last_max = maximum
    draw_progress.last_info_length = wrapped_lines_length

    for _ in range(1 + wrapped_lines_length):
        # this VT100 code causes the console cursor to go up one line
        sys.stdout.write("\x1B[1A")
        # this one clears the line
        sys.stdout.write("\x1B[K")

    # \r at the start of the line causes the new line to overwrite the old
    #   one (in a terminal)
    sys.stdout.write(
        "\r[%s%s] %s\n\r%s\n" % (
            "#" * num_segments,
            " " * (bar_length - num_segments),
            progress_info,
            wrapped_lines
        )
    )
    # if we are at the end, get rid of the abort message
    if progress == 1:
        for _ in range(2 if abortable else 1):
            sys.stdout.write("\x1B[1A")
            sys.stdout.write("\x1B[K")


def wait_for_user(message="continue"):
    """Wait for confirmation by the user (if stdin is interactive)."""
    if sys.stdin.isatty():
        try:
            raw_input("# Press enter to " + message)
        except (KeyboardInterrupt, EOFError):
            # interpret ctrl+d or ctrl+c as an immediate abort
            print  # print one last line for cosmetic purposes
            exit_program()


def times_to_rate_info(times, bunchsize):
    """Convert job completion times to job completion rates.

    bunchsize: jobs are done in batches of this size in parallel

    Returns the mean rate, its error and the minimal rate
    """
    # the first time is always the starting (=base) time
    if len(times) < bunchsize + 1:
        deltas = [float("nan")]
    else:
        # rate = N / time it took to calculate N items, where N is the bunch
        #   size (number of parallel computations)
        # we have to sum over at least that many, because otherwise the
        #   parallel nature will introduce very high rates due to points that
        #   were already finished while others were still written to file
        deltas = [
            (times[i + bunchsize] - times[i]) / bunchsize
            for i in range(len(times) - bunchsize)
        ]

    meandelta = mean(deltas)
    rate = 1 / meandelta
    rateerror = mean([
        abs(1 / meandelta - 1 / (meandelta + p * std(deltas)))
        for p in [-1., 1.]
    ])
    return rate, rateerror, max(rate - rateerror, 1 / max(deltas))


def progress_forecast(current, full, rate):
    """Return estimated completion time based on current progress and rate."""
    if not rate > 0:  # use this instead of <= to also catch nan's
        return "never"
    return time.ctime(time.time() + (full - current) / rate)


def abbreviate(s, width=30):
    """Abbreviate long strings by replacing middle section with '...'."""
    if len(s) <= n:
        return s
    if n < 3:
        return ""
    if n < 5:
        return "..."
    return "%s...%s" % (
        s[0:int(n / 2 - 1)],
        s[int(-(n - n / 2 - 2)):]
    )


def abbreviate_path(path, width=30):
    """Abbreviate too long paths to (preferably unique) shortened ones.

    Uniqueness:
        /folder/file -> /f../file if not other folder starting with "f" exists
    """
    if len(path) <= width:
        return path
    if width < 3:
        return ""

    path = path.replace(os.path.expanduser("~"), "~")
    parts = path.split(os.sep)
    if not parts[0]:
        parts[0] = os.path.abspath(os.sep)

    shortened_parts = [parts[0]]
    for i, name in enumerate(parts[1:-1], start=1):
        parent = os.path.join(*parts[:i])
        siblings = os.listdir(os.path.expanduser(parent))
        for j in range(len(name)):
            shortened = name[:j + 1]
            if not any(sib.startswith(shortened) and
               sib != name for sib in siblings):
                shortened = shortened + ".."
                break
        shortened_parts.append(
            # only save shortened name if actually shorter
            shortened if len(shortened) < len(name) else name
        )
    shortened_parts.append(parts[-1])
    shortened_path = os.path.join(*shortened_parts)

    if len(shortened_path) == width:
        return shortened_path
    elif len(shortened_path) > width:
        for i in range(1, len(shortened_parts)):
            used_parts = ["..."] + shortened_parts[i:]
            shortened_path = os.path.join(*used_parts)
            if len(shortened_path) <= width:
                return shortened_path
        short_filename = abbreviate(parts[-1], width - 4)
        if short_filename in ["...", ""]:
            return "..."
        return ".../" + short_filename
    else:
        best_candidate = ""
        best_usedparts = []
        for k in range(len(shortened_parts) - 1):
            for subset in itertools.combinations(
                    range(len(shortened_parts)), k + 1
                    ):
                used_parts = []
                for i in range(len(shortened_parts)):
                    if i in subset:
                        used_parts.append(shortened_parts[i])
                    else:
                        used_parts.append(parts[i])
                candidate = os.path.join(*used_parts)
                L = len(candidate)
                L0 = len(best_candidate)
                if L <= width and (
                        L > L0 or
                        (L == L0 and max(used_parts) > max(best_usedparts))
                        ):
                    best_candidate = candidate
                    best_usedparts = used_parts
        return best_candidate


def mean(seq):
    """Calculate arithmetic mean of list of numbers."""
    return math.fsum(seq) / len(seq)


def std(seq):
    """Calculate standard deviation of list of numbers."""
    mu = mean(seq)
    return math.sqrt(mean([(x - mu) * (x - mu) for x in seq]))


# ----------------------------- #
# |  [D.04] Iterator Handling | #
# ----------------------------- #

def it_len(iterable):
    """Count elements in an iterator.

    Warning: this consumes the iterator entirely
    """
    c = 0
    for _ in iterable:
        c += 1
    return c


def first(iterable):
    """Return the first element in an iterator."""
    for x in iterable:
        return x


def take(n, iterable):
    """Take n items out of an iterable, starting from the current position."""
    return list(itertools.islice(iterable, n))


def consume(n, iterator):
    """Consume the first n elements in an iterator."""
    # Take the slice starting at position n and consume those before it by
    #   advancing blindly
    next(itertools.islice(iterator, n, n), None)


def roundrobin(*iterables):
    """Return list obtained from interleaving the lists in iterables.

    Example: roundrobin("aaaa", "bbb", "cc") -> a b c a b c a b a a

    Note: this assumes that None is not in any of the iterables
        (satisfied in our use case)
    """
    for xs in itertools.izip_longest(*iterables, fillvalue=None):
        for x in xs:
            if x is not None:
                yield x


def pad_list(iterable, padding, length):
    """Pad list to specified length with padding value."""
    i = -1
    for i, x in enumerate(iterable):
        yield x
    # i + 1 is now length of iterable
    while i + 1 < length:
        yield padding
        i += 1


# ---------------------- #
# |  [D.05] File Input | #
# ---------------------- #

def parspace_file_iterator(file_names):
    """Loop over lines in multiple files and return columns for each line.

    This ignores lines marked as invalid.
    """
    for file_name in file_names:
        with open(file_name) as f:
            for line in f:
                p = map(maybefloat, line.strip().split("\t"))
                if p[0] == ERROR_MARKER:
                    continue
                yield p


def linecount(file_name):
    """Count lines in one or multiple files."""
    if type(file_name) == list:
        return sum(map(linecount, file_name))
    lines = 0
    try:
        if os.path.getsize(file_name) == 0:
            return 0
        # mmap seems to be the fastest way to do this and requires mode "r+"
        with open(file_name, "r+") as f:
            buf = mmap.mmap(f.fileno(), 0)

            rl = buf.readline
            while rl():
                lines += 1
    # ignore any errors regarding missing or unreadable files
    except IOError:
        pass
    except OSError:
        pass
    return lines


def get_columns(columns, row):
    """Extract the columns (given as list of formulas) from a given row."""
    p = []
    for col in columns:
        p.append(
            float(formula_eval(col, file=row))
        )
    return p


def split_cols(cols, ns):
    """Split list into groups with lengths specified by ns.

    This drops anything after sum(ns) columns.
    """
    # if one ns[i] is inf (or NaN), all columns will be funneled into this
    #   one after ones before it are full
    r = [[] for _ in ns]
    i = 0
    for col in cols:
        # use while to skip over ns[i] that are 0
        while len(r[i]) >= ns[i]:
            i += 1

        if i >= len(ns):
            return r

        r[i].append(col)

    return r


def out_repr(x):
    """Convert value to its representation in output files."""
    if isinstance(x, basestring):
        return x.replace("\t", " ").replace("\n", " ")
    return repr(x)


# ----------------------------------------------------------------------------#
# | [D.06] Configuration Input (Parameter Range, Formulas and Preprocessing) |#
# ----------------------------------------------------------------------------#

def find_file(path, base_dir):
    """Find file specified by path.

    It is searched relative to either:
        current folder (1st try)
        base_dir (2nd try)
        or the directory of the running executable (3rd try)
        or the directory of the running executable (if sym link) (3'rd try)
        or relative to the current directory (4th try)

    This always returns an absolute path.
    """
    path = os.path.expanduser(path)
    if os.path.isabs(path):
        return path

    candidate = os.path.abspath(path)
    if os.path.isfile(candidate):
        return candidate

    candidate = os.path.abspath(os.path.join(base_dir, path))
    if os.path.isfile(candidate):
        return candidate

    candidate = os.path.abspath(os.path.join(os.path.dirname(__file__), path))
    if os.path.isfile(candidate):
        return candidate

    if os.path.islink(__file__):
        candidate = os.path.abspath(
            os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                path
            )
        )
        if os.path.isfile(candidate):
            return candidate

    return os.path.abspath(path)


def find_in_path(program):
    """Find executable in PATH environment variable."""
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        if path == ".":
            continue  # ignore "." in PATH since it is not constant
        exe_file = os.path.join(path, program)
        if os.path.exists(exe_file):
            return exe_file

    return None


def find_binary_file(binary_path, base_dir):
    """Find executable file.

    First tries finding it in PATH environment variable, else uses
    find_file function
    """
    # searching through PATH is only necessary for paths that are just names
    if not os.path.isabs(binary_path) and os.path.dirname(binary_path) == "":
        candidate = find_in_path(binary_path)
        if candidate is not None and os.path.isfile(candidate):
            return candidate

    return find_file(binary_path, base_dir)


def PreprocessedConfigFile(file_name, short_name=None, base_dir=None,
                           depth=0, max_depth=10):
    """Return file-like object with @include statements processed."""
    # base_dir is only supposed to be None on the depth = 0 file
    if base_dir is None:
        base_dir = os.path.dirname(find_file(file_name, "."))

    # likewise short_name
    if short_name is None:
        short_name = file_name

    contents = StringIO.StringIO()
    # if not depth:
    #     contents.write("[DEFAULT]\n")
    lineno = 0
    real_locations = []
    has_includes = False
    with open(file_name) as f:
        for line in f:
            lineno += 1
            # format of include lines: "@include" WHITESPACE FILENAME
            if line.startswith("@include"):
                has_includes = True
                if depth + 1 > max_depth:
                    e = ConfigParser.ParsingError(file_name)
                    e.append(
                        lineno,
                        "maximal include depth exceeded: %i" % (depth + 1)
                    )
                    raise e

                # strip the trailing "\n" and split into the one "@include"
                #   and the rest of the line with arbitrary whitespace after
                #   the "@include"
                other_filename = line.rstrip("\n").split(None, 1)[-1]
                full_name = find_file(other_filename, base_dir)

                other_file = PreprocessedConfigFile(
                    full_name,
                    short_name=other_filename,
                    base_dir=base_dir,
                    depth=depth + 1
                )
                real_locations.extend(other_file.real_locations)
                for line in other_file:
                    contents.write(line)

                # make sure the next line in the including file starts normally
                if not other_file.getvalue().endswith("\n"):
                    contents.write("\n")

            else:
                real_locations.append((lineno, short_name))
                contents.write(line)

    # context support: translate ConfigParser error if appropriate
    setattr(contents, "__enter__", lambda contents=contents: contents)

    def exit_preproc_context(exc_type, exc_value, _tb, contents=contents):
        contents.close()
        if exc_type == ConfigParser.ParsingError:
            raise contents.translate_error(exc_value)
    setattr(contents, "__exit__", exit_preproc_context)

    # setup the real locations info
    # error translation is only necessary if includes happened
    setattr(contents, "real_locations", real_locations)
    if has_includes:
        def translate_error(e, self=contents):
            # ParsingError.errors has format: [(lineno, text), ...]
            e.message = e.message.splitlines()[0]
            for error_data in e.errors:
                e.message += '\n\t[line %i in %s]: %s' % (
                    self.real_locations[error_data[0] - 1][0],
                    self.real_locations[error_data[0] - 1][1],
                    error_data[1]
                )

            return e
        contents.translate_error = translate_error
    else:
        contents.translate_error = lambda e: e

    contents.seek(0)
    return contents


# namedtuple to differentiate between different range types and save the
#   associated range parameters
# NOTE: these classes have to be picklable, so we have to define them at
#   top level
ParameterRange_Interval = namedtuple(
    "ParameterRange_Interval", ["start", "end"]
)
ParameterRange_NormalVariate = namedtuple(
    "ParameterRange_NormalVariate", ["mu", "sigma"]
)


class ParameterRange(object):

    """Class used to convey information on parameter ranges.

    They can be given in the following form:
        * discrete/finite ranges: "a, b, c" and so on, or "a, b, ..., c"
            (meaning the ellipsis is substituted like "a, b, a+2(b-a), ..., c"
            these can also be concatenated, e.g. "1, 2, ..., 10, 20, ..., 100"
                -> (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70,
                    80, 90, 100)
            alternative syntax: "finite(a, b, ... c)" etc.
        * continuous ranges:
            * on an interval:
                "interval(a, b) with [count=<n>,] [distribution=<log|linear>]
                if count is given, the range will automatically be converted
                    to a finite range with spacing according to distribution
                    (default:linear)
            * following normal distribution:
                "normalvariate(mu, sigma) with [count=<n>]"
                similarly to interval, specifying count makes the range
                    finite while preserving its distribution as sample

        options and type names are case-insensitive
    This class uses decimal internally, but all exposed values are float.
    """

    def __init__(self, s):
        """Initialize parameter range with definition."""
        if not isinstance(s, basestring):
            raise ValueError(
                "invalid parameter range definition type: " + str(type(s))
            )
        if not s:
            raise ValueError("empty parameter range")
        parts = s.split("with")
        if len(parts) > 2:
            raise ValueError("too many option parts: " + s)

        range_def, options = map(str.strip, (parts + [""])[0:2])
        if not range_def:
            raise ValueError("empty definition: " + s)

        try:
            self._options = dict()
            for k_v in options.split(","):
                if k_v:
                    k, v = (map(str.strip, k_v.split("=")) + [""])[0:2]
                    # Convert with .lower() for case-insensitivity
                    self._options[k.lower()] = maybedecimal(v)
        except Exception as e:
            raise ValueError("error parameter range options: " + str(e))

        if range_def[0] in string.ascii_letters:
            iterator = iter(range_def)
            # Convert with .lower() for case-insensitivity
            type_name = ''.join(
                itertools.takewhile(lambda x: x != "(", iterator)
            ).lower()
            definition = ''.join(
                itertools.takewhile(lambda x: x != ")", iterator)
            )
            garbage = ''.join(iterator)
            if garbage:
                raise ValueError(
                    "unexpected part in range definition: " + garbage
                )
        else:
            type_name = "finite"
            definition = range_def

        if type_name == "interval":
            self._type = "interval"
            if "distribution" not in self._options:
                self._options["distribution"] = "linear"

            parts = definition.split(",")
            if len(parts) != 2:
                raise ValueError("invalid interval definition: " + definition)
            try:
                start, end = map(Decimal, parts)
                if start == end:
                    raise ValueError("zero interval length")
                if start > end:
                    raise ValueError("negative interval length")
                self._definition = "interval(%s, %s)" % (start, end)
            except (ValueError, decimal.InvalidOperation) as e:
                raise ValueError("invalid interval bounds: " + str(e))

            dist = self._options.get("distribution", "linear")
            if dist not in ["linear", "log"]:
                raise ValueError("unsupported distribution: " + str(dist))

            if dist == "log" and start <= 0 or end <= 0:
                raise ValueError(
                    "invalid bounds for log distribution: %s, %s" % (
                        start, end
                    )
                )

            if "count" not in self._options:
                self._data = ParameterRange_Interval(
                    start=float(start), end=float(end)
                )
            else:
                self._type = "finite"
                count = int(self._options["count"])
                if count < 2:
                    raise ValueError("point count too small: " + str(count))

                if dist == "linear":
                    self._data = []
                    delta = (end - start) / (count - 1)

                    for i in xrange(count - 1):
                        self._data.append(start + i * delta)
                    # add end by hand instead of calculating it due to
                    #   rounding errors
                    self._data.append(end)

                elif dist == "log":
                    self._data = []

                    delta = (end / start) ** (1 / Decimal(count - 1))
                    for i in xrange(count - 1):
                        self._data.append(start * (delta ** i))
                    # some numerical errors always creep in, so add the end
                    #   by hand and not by calculation
                    self._data.append(end)

        elif type_name == "normalvariate":
            self._type = "normalvariate"
            parts = definition.split(",")
            if len(parts) != 2:
                raise ValueError("invalid interval definition: " + definition)
            try:
                mu, sigma = map(Decimal, parts)
                if sigma <= 0:
                    raise ValueError("non-positive standard deviation")
                self._data = ParameterRange_NormalVariate(
                    mu=float(mu), sigma=float(sigma)
                )
                self._definition = "normalvariate(%s, %s)" % (mu, sigma)
            except (ValueError, decimal.InvalidOperation) as e:
                raise ValueError("invalid normalvariate parameters: " + str(e))

            if "distribution" in self._options:
                raise ValueError(
                    "'distribution' option not supported for normalvariate"
                )

            if "count" in self._options:
                self._type = "finite"
                count = int(self._options["count"])
                if count < 1:
                    raise ValueError("point count too small: " + str(count))
                self._data = []
                delta = 1.0 / (count + 1)
                for i in xrange(count):
                    self._data.append(
                        float(mu) +
                        float(sigma) * inversenormalcdf(delta * (i + 1))
                    )

        elif type_name == "finite":
            self._type = "finite"
            if "distribution" in self._options:
                raise ValueError(
                    "'distribution' option not supported for finite ranges"
                )
            if "count" in self._options:
                raise ValueError(
                    "'count' option not supported for finite ranges"
                )
            parts = map(maybedecimal, map(str.strip, definition.split(",")))
            self._data = []
            self._definition = ', '.join(map(str, parts))
            for i, x in enumerate(parts):
                if x in ["..", "..."]:
                    if i - 2 < 0 or i + 1 >= len(parts) or \
                       not is_number(parts[i - 2]) or \
                       not is_number(parts[i - 1]) or \
                       not is_number(parts[i + 1]):
                        # any of the parts is not there or not a number
                        raise ValueError(
                            "ellipsis must be preceded by two valid numbers "
                            "and succeeded by one"
                        )
                    # range of form a, b, ... ,c
                    a, b, c = (
                        Decimal(parts[i - 2]),
                        Decimal(parts[i - 1]),
                        Decimal(parts[i + 1])
                    )
                    if a == b:
                        raise ValueError(
                            "range with step 0: %s == %s" % (a, b)
                        )
                    if not (a < b < c or a > b > c):
                        raise ValueError(
                            "unordered partial range: does not fulfill "
                            "%s < %s < %s" % (a, b, c)
                        )

                    delta = b - a
                    x = b + delta
                    while x < c:
                        self._data.append(float(x))
                        x += delta
                else:
                    try:
                        x = Decimal(x)
                    except (ValueError, decimal.InvalidOperation) as e:
                        raise ValueError("invalid range part: " + str(e))

                    self._data.append(float(x))
        else:
            raise ValueError("unsupported range type: " + type_name)

        if self._type == "finite":
            self._data = map(float, self._data)
            self._strdata = map(repr, self._data)

        if "distribution" in self._options:
            self._options["distribution"] = str(self._options["distribution"])

        if "mcmc_stepsize" in self._options:
            self._options["mcmc_stepsize"] = float(
                self._options["mcmc_stepsize"]
            )

    def __str__(self):
        """Convert parameter range back to definition string."""
        if self._options:
            return self._definition + " with " + ', '.join(
                sorted(
                    k + "=" + str(v)
                    for k, v in self._options.iteritems()
                )
            )
        return self._definition

    def __repr__(self):
        """Return evaluable object definition."""
        return "ParameterRange(%r)" % str(self)

    @property
    def is_finite(self):
        """Return whether parameter range has only finite amount of values."""
        return isinstance(self._data, list)

    @property
    def values(self):
        """Return tuple of possible parameter values."""
        if self.is_finite:
            return tuple(self._data)
        else:
            return self._data

    @property
    def options(self):
        """Return the options used in parameter range definition."""
        return dict(self._options)

    def index(self, value):
        """Return the index of value in our list of values."""
        if isinstance(self._data, list):
            return self._strdata.index(repr(value))
        raise ValueError(
            "parameter range type '%s' is not iterable" % self._type
        )

    def truncate(self, value):
        """Truncate value to one of the possible values.

        This will be the one with the same repr representation.
        """
        if isinstance(self._data, list):
            return self._data[self._strdata.index(repr(value))]
        raise ValueError(
            "parameter range type '%s' is not iterable" % self._type
        )

    def pull_inside(self, value):
        """Return number in parameter range that is nearest to given value."""
        if isinstance(self._data, list):
            return min(self._data, key=lambda x: abs(x - value))
        elif isinstance(self._data, ParameterRange_NormalVariate):
            return value
        elif isinstance(self._data, ParameterRange_Interval):
            if self._data.start < self._data.end:
                return max(self._data.start, min(value, self._data.end))
            else:
                return max(self._data.end, min(value, self._data.start))

    def __nonzero__(self):
        """Return True."""
        # non-zero length parameter ranges are filtered out by __init__
        # thus parameter ranges are always considered non-zero
        return True

    def __len__(self):
        """Return number of possible values (if finite)."""
        if isinstance(self._data, list):
            return len(self._data)
        raise ValueError(
            "parameter range type '%s' has no length" % self._type
        )

    def __contains__(self, value):
        """Return whether value is member of parameter range."""
        if isinstance(self._data, list):
            # values are often converted to string and back to float
            # they should still be "inside" the values list afterwards
            return repr(value) in self._strdata
        elif isinstance(self._data, ParameterRange_NormalVariate):
            return True
        elif isinstance(self._data, ParameterRange_Interval):
            if self._data.start < self._data.end:
                return self._data.start < value < self._data.end
            else:
                return self._data.end < value < self._data.start

    def __iter__(self):
        """Return iterator over parameter range values."""
        if isinstance(self._data, list):
            return iter(self._data)
        raise ValueError("parameter range type '%s' not iterable" % self._type)


# namedtuple like class that does not care about the number of names it gets
#   the names label the first len(names) entries, the rest can only be
#   accessed by index
# NOTE: member access is only about 30% slower than in namedtuple
def flexible_record(names):
    """Return namedtuple-like object that does not care about names."""
    class _flexible_record(tuple):
        pass

    _flexible_record.__names__ = []
    for index, name in enumerate(names):
        _flexible_record.__names__.append(name)
        # add property to new record class that references indexed element
        #   by name
        setattr(
            _flexible_record,
            name,
            property(lambda self, index=index: self[index])
        )

    return _flexible_record


def check_names(names):
    """Check whether the given list yields valid names.

    That means no empty names, no duplicates, no keywords,
    no invalid characters.
    """
    names_set = set()
    identifier_chars = string.ascii_letters + string.digits + "_"
    for n in names:
        if not n:
            raise ValueError("empty name used")
        if n in names_set:
            raise ValueError("name used multiple times: %s" % n)
        if keyword.iskeyword(n):
            raise ValueError("name is keyword: %s" % n)
        if n.startswith("_") or not all(c in identifier_chars for c in n):
            raise ValueError("name contains invalid characters: %s " % n)
        names_set.add(n)
    return


class FormulaError(Exception):

    """Exception class raised in formulas by formula_eval."""

    def __str__(self):
        """Return concatenated error arguments.

        First argument is the error message.
        Second is parameter point that generated the error.
        """
        return '\n'.join(self.args)


def formula_eval(function,
                 pars=None, vars=None,
                 data=None, file=None,
                 **kwargs):
    """Compile and evaluate user supplied Python expression formulas."""
    # Has to be called as in
    #   formula_eval(None, par_names, var_names, data_names,
    #       file_column_names)
    #   to set names
    if not hasattr(formula_eval, "cache"):
        formula_eval.cache = {}

    cache = formula_eval.cache
    if function is None:
        # cache[None] is used to store the namedtuple converter
        # arguments now specify the names for the different values
        # NOTE: names have to be valid identifiers:
        #   letters, digits, and underscores
        #   not starting with underscore
        #   no Python keyword
        cache[None] = (
            namedtuple("par_record", pars),
            namedtuple("var_record", vars),
            # data values, file columns can also be wholly unnamed if the
            #   user does not specify names
            # -> no need to use flexible record then
            flexible_record(data) if data else (lambda x: x),
            flexible_record(file) if file else (lambda x: x)
        )
        # namedtuple will raise ValueError's if one of the names is not valid
        #   should never happen due to check_names called before this
        return

    par_record, var_record, data_record, file_record = cache[None]

    # key is code + used keyword args, such that the compiled code is uniquely
    #   determined by the call
    key = "".join([function, ":"] + kwargs.keys())
    if key in cache:
        cached_function = cache[key]
    else:
        # replace new line chars, so the return is definitely the only
        #   executed statement
        code = "def function(pars, vars, data, file%s): return (%s)" % (
            ', ' + ', '.join(kwargs.keys()) if kwargs else '',
            function.replace("\n", " ")
        )
        namespace = math_namespace()
        namespace.update(formula_eval.helper_modules)
        # this inherits future division implicitly
        exec code in namespace
        cached_function = cache[key] = namespace["function"]
    try:
        # only pars, vars, data and file items have names everything else is
        #   only by index or not iterable
        return cached_function(
            pars=par_record(*pars) if pars else None,
            vars=var_record(*vars) if vars else None,
            data=data_record(data) if data else None,
            file=file_record(file) if file else None,
            **kwargs
        )
    except Exception as e:
        raise FormulaError(
            "Error in formula '%s': [%s] %s" % (
                function.replace("\n", " "),
                type(e).__name__,
                str(e)
            ),
            "For pars=%s, vars=%s, data=%s, file=%s%s" % (
                pars, vars, data, file,
                (
                    ", " + ", ".join(
                        k + "=" + repr(v) for k, v in kwargs.iteritems()
                    )
                ) if kwargs else ""
            )
        )


def remember(*args, **kwargs):
    """Return or save named buffered calculation results."""
    if not args:
        if len(kwargs) != 1:
            raise ValueError(
                "remember called with more or less than one parameter"
            )
        remember.memory.update(kwargs)
        return kwargs.values()[0]
    if args[0] is None:
        remember.memory = {}
        return
    return remember.memory[args[0]]


def math_namespace():
    """Return dictionary containing math functions for formulas."""
    return {
        "math": math,
        "pi": math.pi,

        "cos": math.cos, "sin": math.sin, "tan": math.tan,
        "acos": math.acos, "asin": math.asin,
        "atan": math.atan, "atan2": math.atan2,

        "exp": math.exp, "log": math.log, "log10": math.log10,
        "sqrt": math.sqrt,
        "acosh": math.acosh, "asinh": math.asinh, "atanh": math.atanh,
        "cosh": math.cosh, "sinh": math.sinh, "tanh": math.tanh,
        "remember": remember
    }


# ---------------------- #
# |  [D.07] Randomness | #
# ---------------------- #

def scatter_iterator(count, par_ranges):
    """Yield a given number of random points in the given parameter ranges."""
    for _ in xrange(count):
        yield make_random_point(par_ranges)


def make_random_point(par_ranges):
    """Return random point with the distribution of the parameter ranges."""
    point = []
    # increase speed by minimizing the amount of look-ups
    append = point.append
    randomuniform = random.uniform
    randomnormalvariate = random.normalvariate
    randomchoice = random.choice

    for r in par_ranges:
        values = r.values
        # if parameter range is continuous, return new point with appropriate
        #   distribution
        if not r.is_finite:
            if isinstance(values, ParameterRange_Interval):
                dist = r.options.get("distribution", "linear")
                if dist == "linear":
                    append(randomuniform(values.start, values.end))
                elif dist == "log":
                    append(
                        math.exp(
                            randomuniform(
                                math.log(values.start),
                                math.log(values.end)
                            )
                        )
                    )
                else:
                    raise ValueError("unsupported distribution type: " + dist)
            elif isinstance(values, ParameterRange_NormalVariate):
                append(randomnormalvariate(values.mu, values.sigma))
            else:
                raise ValueError(
                    "unsupported distribution type: " + str(type(values))
                )
        # otherwise jump randomly in discrete list with appropriate step size
        else:
            append(randomchoice(values))

    return point


# ---------------------------------------------- #
# |  [D.08] Markov Chain Monte Carlo Functions | #
# ---------------------------------------------- #

class ProbabilityError(ValueError):

    """Custom exception for the Metropolis-Hastings algorithm.

    Signals that a candidate point has vanishing likelihood
    """

    pass


def allowed_by_prior(x, par_ranges):
    """Check whether given point has non-zero prior likelihood."""
    for i, r in enumerate(par_ranges):
        if x[i] not in r:
            return False
    return True


def prior_likelihood(x, par_ranges):
    """Calculate prior likelihood for given point and parameter ranges."""
    p = 1
    for i, r in enumerate(par_ranges):
        # contained check also works for continuous ranges
        if x[i] not in r:
            p = 0
            break

        values = r.values
        if isinstance(values, ParameterRange_Interval):
            dist = r.options["distribution"]
            if dist == "linear":
                pass  # p *= 1
            elif dist == "log":
                # dp = C dlogx = C 1/x * dx
                p *= 1 / abs(x[i])
            else:
                raise ValueError("unsupported distribution type: " + str(dist))
        elif isinstance(values, ParameterRange_NormalVariate):
            p *= math.exp(-0.5 * ((x[i] - values.mu) / values.sigma) ** 2)
        elif isinstance(values, tuple):
            # finite -> Laplace probability (all equal)
            pass  # p *= 1
        else:
            raise ValueError(
                "unsupported distribution type: " + str(type(values))
            )

    return p


def mcmc_step(theta0, par_ranges):
    """Make one iteration "step".

    That means: return new point theta, obtained from proposal probability
        density q(theta, theta0) around theta0, together with the value
        of q(theta, theta0) and q(theta0, theta)
    """
    theta, q, q0 = [], 1, 1
    append = theta.append
    randomnormalvariate = random.normalvariate
    for i, r in enumerate(par_ranges):
        if "mcmc_stepsize" in r.options:
            # standard Gaussian stepsize -> symmetric
            # is always canceled in Hastings test ratio
            # q *= 1
            # q0 *= 1
            if r.is_finite:
                # if there are infinite indices, q(i1,i0) is symmetric
                i0 = r.index(theta0[i])
                i1 = int(round(
                    randomnormalvariate(i0, r.options["mcmc_stepsize"])
                ))

                # however, if we go out of bounds on the new index, the
                #   implicit prior on the interval [0, len(r.values)] takes
                #   effect and rules out the point
                # raising an exception to the MH algorithm correctly handles
                #   this
                if i1 < 0 or len(r) <= i1:
                    raise ProbabilityError("out of finite bounds")
                append(r.values[i1])
            else:
                append(
                    randomnormalvariate(theta0[i], r.options["mcmc_stepsize"])
                )
        else:
            # take from prior density -> not-symmetric but well known
            x = make_random_point([r])[0]
            q *= prior_likelihood([x], [r])
            q0 *= prior_likelihood([theta0[i]], [r])
            append(x)

    return (theta, q, q0)


def make_chain(data):
    """Standard (Metropolis-)Hastings algorithm to calculate a MCMC sample.

    Algorithm correspondence:
        p(theta) = likelihood(theta) * prior(theta) * has_errors(theta)
        q(theta, theta0) = separable into densities for one coordinate each:
            if mcmc_stepsize defined:
                normalvariate around theta0 = f( (theta-theta0)_i^2 )
                -> symmetric
            else:
                q(theta, theta0) = prior(theta) -> not symmetric
                q(theta, theta0) and q(theta0, theta) are returned by mcmc_step
    """
    (
        output_file, rejected_file, x0, current_length, final_length,
        par_ranges, current_iteration, likelihood_code, heating_function_code,
        status_file_path, debug_mode, random_seed
    ) = data

    def likelihood(x, code=likelihood_code):
        return formula_eval(code, pars=x[0], vars=x[1], data=x[2])

    if heating_function_code:
        def heating_function(iteration, dlogL, logL0):
            return formula_eval(
                heating_function_code,
                iteration=iteration,
                dlogL=dlogL,
                logL0=logL0
            )

    else:
        heating_function = None

    L0 = likelihood(x0) * prior_likelihood(x0[0], par_ranges)

    # IMPORTANT: re-seed PRNG so all chains are guaranteed to find distinct
    #   points - THIS HAS TO BE DIFFERENT FOR EACH CHAIN
    random.seed(random_seed)

    staycount = 0
    # we keep track of all errors and reasons why points get discarded
    #   this resets when new point is found and is only ever shown in expert
    #   mode
    reasons = {"prior": 0, "likelihood": 0, "errors": 0, "chance": 0}
    try:
        while current_length < final_length:
            with open(status_file_path, "w") as status_file:
                pickle.dump([
                    current_length,
                    current_iteration,
                    reasons
                ], status_file)

            current_iteration += 1
            alpha = 0
            try:
                # x0[0] because x is (pars, vars, data)
                x1, q1, q0 = mcmc_step(x0[0], par_ranges)
                if not allowed_by_prior(x1, par_ranges):
                    x1 = Exception(x1, "excluded by prior = 0")
                    reasons["prior"] += 1
                    staycount += 1
                    raise ProbabilityError("point excluded by prior = 0")
                x1 = process_item(x1)
                if isinstance(x1, Exception):
                    reasons["errors"] += 1
                    staycount += 1
                    raise ProbabilityError("point excluded by processor error")

                L1 = likelihood(x1) * prior_likelihood(x1[0], par_ranges)

                if L1 <= 0:
                    reasons["likelihood"] += 1
                    staycount += 1
                    raise ProbabilityError("point excluded by likelihood = 0")

            except ProbabilityError as e:
                # all other errors are not recoverable and break the chain
                #   immediately
                L1 = 0

            else:  # executed if try block was finished without exception
                # Hastings testing ratio
                #        p(theta)  q(theta0, theta)
                # min(1, -------------------------- )
                #        p(theta0) q(theta, theta0)
           
                alpha = min(
                    1,
                    (
                        L1 / L0 if heating_function is None else math.exp(
                            heating_function(
                                current_iteration,
                                math.log(L1 / L0),
                                math.log(L0)
                            )
                        )
                    ) * q0 / q1
                )

# Edit added by Ciara - missing step where all points for which L1>=1 should 
# be immiediately accepted, added 'or' statement to fix this
            
            if random.random() < alpha or L1 >= 1:
                with open(output_file, "a") as f:
                    f.write(
                        "\t".join(
                            map(out_repr, list_sum(x0) + [L0, staycount])
                        )
                    )
                    f.write("\n")
                x0 = x1
                L0 = L1
                staycount = 1
                current_length += 1
                reasons = {
                    "prior": 0, "likelihood": 0, "errors": 0, "chance": 0
                }
            else:
# Warning/WARNING/
#               if not isinstance(x1, Exception):
#                   with open(rejected_file, "a") as f:
#                       f.write("\t".join(map(out_repr, list_sum(x1) + [L1])))
#                       f.write("\n")

                staycount += 1
                if alpha > 0:
                    reasons["chance"] += 1

        # write the last point to file
        with open(output_file, "a") as f:
            f.write("\t".join(map(out_repr, list_sum(x0) + [L0, staycount])))
            f.write("\n")
        # and write the 100% status, too
        with open(status_file_path, "w") as status_file:
            pickle.dump([
                current_length,
                current_iteration,
                reasons
            ], status_file)

    except Exception as e:
        if debug_mode:
            traceback.print_exc()
        raise Exception("Chain ended by error: " + str(e))

    return "Chain done in file %r" % file
    # nobody will ever read this return value anyway


# ----------------------------------- #
# |  [D.09] Explorer Mode Functions | #
# ----------------------------------- #

def index2point(k, par_ranges):
    """Convert index into parameter space grid into coordinates.

    The parameter space grid is counted the following way:
        Let X_i = { x_ij } be the set of values for parameter i
        Then the index k of (x_{1,n_1}, x_{2,n_2}, x_{3,n_3}, ...) is
            k = n1 + |X_1| n2 + |X_1| |X_2| n3 + ...

        Then
            n_1 = k mod |X_1|
            n_2 = (k - n_1) / |X_1| mod |X_2|
            n_3 = (k - n_1 - |X_1| n_2) / (|X_2| |X_1|) mod |X_3|
            ...
    """
    # Variables:
    # n_i = ki[i - 1] (arrays start at 0)
    # k_i = (k - km) / l mod |X_i|
    #   (with km and l defined such that it matches last equations in doc str

    km = 0
    l = 1
    ki = [0] * len(par_ranges)
    for i, r in enumerate(par_ranges):
        ki[i] = int((k - km) / l % len(r))
        km += ki[i] * l
        l *= len(r)
    return tuple([r[ki[i]] for i, r in enumerate(par_ranges)])


def point2index(x, par_ranges):
    """Inverse function of index2point."""
    # Note that range(0) == [] and prod([]) == 1
    return sum(
        par_ranges_i.index(x[i])
        *
        prod([len(par_ranges[j]) for j in range(i)])
        for i, par_ranges_i in enumerate(par_ranges)
    )


def neighbors(index, par_ranges, stepsize=1):
    """Return neighbors of point.

    Neighbors defined as in:
      N
    N X N
      N
    This scales as 2*d for d dimensions
    """
    result = []
    rs = [range(len(r)) for r in par_ranges]
    multiindex = index2point(index, rs)
    stepsize_ = int(max(math.ceil(stepsize), 1))
    for i, r in enumerate(par_ranges):
        neighbor_multiindex = list(multiindex)
        neighbor_multiindex[i] = min(multiindex[i] + stepsize_, len(r) - 1)
        if neighbor_multiindex[i] != multiindex[i]:
            result.append(tuple(neighbor_multiindex))
        neighbor_multiindex[i] = max(multiindex[i] - stepsize_, 0)
        if neighbor_multiindex[i] != multiindex[i]:
            result.append(tuple(neighbor_multiindex))
    return [point2index(n, rs) for n in result]


def extended_neighbors(index, par_ranges, stepsize=1):
    """Return extended neighbors of given point.

    Neighbors defined as in
    N N N
    N X N
    N N N

    This scales as 3^d-1 for d dimensions
    """
    rs = [range(len(r)) for r in par_ranges]
    multiindex = index2point(index, rs)
    stepsize_ = int(max(math.ceil(stepsize), 1))

    dim = len(multiindex)
    result = [
        tuple(
            [
                min(max(multiindex[i] + delta[i], 0), len(par_ranges[i]) - 1)
                for i in range(dim)
            ]
        )
        for delta in itertools.product([-stepsize_, 0, stepsize_], repeat=dim)
    ]
    result.remove(multiindex)
    return [point2index(n, rs) for n in result]


def is_boundary(index, nodes, par_ranges):
    """Check if a point lies on the boundary of the set of points in nodes.

    A point lies on the boundary if and only if:
        * it has been calculated and was valid
        * it is "likely" enough
        * one of its neighbors has not been calculated yet
    """
    index_data = nodes[index]
    # nodes[index][1] is 'maybe_boundary' field
    #   it is initially True and is set to False (only here!) when all
    #   neighbors have been calculated
    if index_data is None or not index_data[1] or \
       not index_data[0] >= min_likelihood:
        # use "not a >= b" instead of "a < b" to catch "nan"s
        return False

    for n in neighbors(index, par_ranges):
        if n not in nodes:
            return True

    index_data[1] = False
    return False


def interpolate_indices(i1, i2, par_ranges):
    """Return the points in parameter space connecting two points."""
    if i1 == i2:
        return []
    value_indices = [range(len(r)) for r in par_ranges]
    # first convert index to multiindex where each component is the index
    #   into the set of parameter values
    v1 = index2point(i1, value_indices)
    v2 = index2point(i2, value_indices)
    # calculate the interpolation direction...
    delta = [v2[i] - v1[i] for i in range(len(par_ranges))]
    # ...and the coordinate direction in which we have the most steps
    # this will parametrize the line on parameter space grid
    longest_dir = max(range(len(par_ranges)), key=lambda i: abs(delta[i]))

    return [
        point2index(p, value_indices) for p in [
            [
                round(v1[i] + (delta[i] * 1.0 * n) / abs(delta[longest_dir]))
                for i in range(len(par_ranges))
            ]
            # take all n from 0 to N -> range(N+1)
            for n in range(abs(delta[longest_dir]) + 1)
        ]
    ]


def extrapolate_indices(i1, i2, par_ranges):
    """Extrapolate along the line defined by two points.

    That means for point x1 indexed by i1 and x2 indexed by i2
        Find the new point x3 = x1 + 2 (x2 - x1) and all points between
        x2 and x3
    """
    value_indices = [range(len(r)) for r in par_ranges]
    v1 = index2point(i1, value_indices)
    v2 = index2point(i2, value_indices)

    v3 = [
        max(0, min(len(par_ranges[i]) - 1, v1[i] + 2 * (v2[i] - v1[i])))
        for i in range(len(v1))
    ]

    return interpolate_indices(i2, point2index(v3, value_indices), par_ranges)


# ---------------------------- #
# |  [D.10] Point Processing | #
# ---------------------------- #

def process_item(point, test_mode=False, print_error_trace=False):
    """Process the point given by the values in <point>.

    Result is either a tuple (pars, vars, data) or an Exception (as an object)
    The exception error message is the first caught exception and most likely
        the first bound that is not satisfied.
    """
    # we get an already initialized point_processor from the main process
    #   through os.fork()
    global point_processor
    # everything else is transferred via attributes of process_item

    (
        pars_with_names, vars_with_names, par_and_vars_with_names
    ) = process_item.with_names
    point = tuple(point)
    remember(None)

    result = Exception(point, "Unknown Error")
    # make a temporary directory in which we can run our item processor
    #   without any disturbances
    # keep the dir around if requested (i.e. in test mode)
    with TemporaryDirectory(keep_dir=test_mode) as tmpdir:
        if test_mode:
            print "# Working directory:", tmpdir
        with pushd(tmpdir):
            try:
                point_vars = []
                for formula in process_item.var_formulas:
                    point_vars.append(formula_eval(formula, pars=point))

                dressed_pars = pars_with_names(*point)
                dressed_vars = vars_with_names(*point_vars)
                dressed_both = par_and_vars_with_names(
                    *(tuple(point) + tuple(point_vars))
                )

                processed_template_file = None
                if process_item.template:
                    processed_template_file = os.path.join(tmpdir, "template")
                    with open(processed_template_file, "w") as \
                            processed_template:
                        processed_template.write(
                            eval(
                                process_item.template,
                                {
                                    "pars": dressed_pars,
                                    "vars": dressed_vars,
                                    "values": dressed_both
                                },
                                dressed_both._asdict()
                            )
                        )

                # if the processor takes more than one argument, it also
                #   expects pars and vars as namedtuples
                if point_processor and \
                   point_processor.main.func_code.co_argcount == 1:
                    processed_point = point_processor.main(
                        processed_template_file
                    )
                elif point_processor:
                    processed_point = point_processor.main(
                        processed_template_file,
                        dressed_pars,
                        dressed_vars
                    )
                else:
                    processed_point = []

                check_bounds(
                    process_item.bounds,
                    point, point_vars, processed_point
                )

# - Warning - made a change here
                return (list(point), list(point_vars), list(processed_point))
#                return (list(processed_point))

            except Exception as e:
                if print_error_trace:
                    traceback.print_exc()
                # For FormulaError, only the first arg is relevant, since the
                #   second one is the point for which it was evaluated and we
                #   already have that at hand
                if isinstance(e, FormulaError):
                    e = Exception(e.args[0])

                # if the exception has an output attribute, it is probably
                #   and CalledProcessError, which is interesting for debugging
                if test_mode and hasattr(e, "output"):
                    print "# Error output:\n", e.output

                # make sure the error message stays in one column and line by
                #   replacing \t and \n
                msg = str(e).replace("\t", " ").replace("\n", " ")
                return Exception(point, msg)

    return result


def process_annotated_item(point_and_info):
    """Process point and keep track of info attached to it."""
    (point, info) = point_and_info
    data = process_item(point)
    return (data, info)


def process_batch(batch):
    """Process multiple points at once."""
    if not batch:
        return []
    # Sequence means it has __getitem__ and __len__, so batch[0][0] works
    if not isinstance(batch[0], Sequence):
        raise ValueError("Invalid batch")

    # if batch has structure [(Sequence, something), ...], each list item
    #   must have form (parameter_values, annotation)
    if isinstance(batch[0][0], Sequence):
        result = process_batch.pool.imap(
            process_annotated_item, batch
        )
    else:
        result = process_batch.pool.imap(process_item, batch)

    return result


class Bunch:

    """Object with a bunch of named attributes."""

    def __init__(self, **kwds):
        """Construct the object."""
        self.__dict__.update(kwds)

    __str__ = __repr__ = lambda self: 'Bunch(**%s)' % repr(dict(self.__dict__))


def process_batch_as_worker(batch):
    """Process multiple points and keep track of the progress."""
    if not batch:
        raise StopIteration()
    batch_length = len(batch)
    start_time = time.time()
    batch_id = hash(start_time)

    process_batch.running_batches[batch_id] = Bunch(
        current=0,
        full_length=batch_length,
        start_time_str=time.strftime("%d.%m.%Y %H:%M:%S"),
        start_time=start_time
    )
    log_msg("# Starting batch (id %s) of %i points" % (batch_id, batch_length))
    try:
        for i, x in enumerate(process_batch(batch)):
            process_batch.running_batches[batch_id].current = i + 1
            yield x
    finally:
        # the last iteration does not return after its yield, so use a
        #   try/finally block
        log_msg("# Finished batch (id %s)" % batch_id)


def check_bounds(bounds, pars, vars, data):
    """Check all bounds and raise ValueError if one is violated."""
    for bound in bounds:
        allowed = formula_eval(
            bound,
            pars=pars,
            vars=vars,
            data=data
        )
        if not allowed:
            raise ValueError("point violates: %s" % bound)


def apply_symmetry(rule, point, par_names):
    """Apply symmetry to parameter point."""
    transformed_point = list(point)
    replaced_values = formula_eval(
        rule,
        pars=transformed_point, **dict(zip(par_names, par_names))
    ).iteritems()
    for k, v in replaced_values:
        if isinstance(k, int):
            transformed_point[k] = v
        elif isinstance(k, basestring):
            transformed_point[par_names.index(k)] = v

    return tuple(transformed_point)


def processing_concurrency():
    """Return local processing concurrency.

    This is used for weighting of worker processes.
    """
    return processing_concurrency.value


# -------------------------------- #
# |  [D.11] Program Flow Control | #
# -------------------------------- #

def init_processing(
        par_names, var_names, data_names, file_column_names,
        var_formulas, bounds, template_file, config_dir, helper_modules
):
    """Initialize formula_eval, process_item and the template system."""
    # example: helper_modules = /path/to/module_a.py:/other/path/to/module_b.py
    # -> formulas have module_a and module_b objects available
    formula_eval.helper_modules = dict()
    for module in helper_modules:
        try:
            module_path = find_file(module, config_dir)
            if os.path.isfile(module_path):
                module_name = os.path.splitext(
                    os.path.basename(module_path)
                )[0]
                module_obj = imp.load_source(module_name, module_path)
            else:
                # if it's not a file, maybe it's something Python finds itself
                # only take base module name to make usage of e.g.
                #   scipy.optimize possible when specifying it directly
                module_name = module.split(".")[0]
                # also if you import "a.b", Python return 'a' instead
                module_obj = __import__(module)
            formula_eval.helper_modules[module_name] = module_obj
        except Exception as e:
            raise ImportError("could not import '%s': %s" % (module, e))

    # this gives formula_eval the names and does no actual evaluation
    formula_eval(
        None,
        pars=par_names, vars=var_names, data=data_names,
        file=file_column_names
    )

    # give process_item its static information:
    # this makes it possible to access pars and vars members via pars.NAME
    #   and vars.NAME
    process_item.with_names = (
        namedtuple("pars_with_names", par_names),
        namedtuple("vars_with_names", var_names),
        namedtuple(
            "par_and_vars_with_names", par_names + var_names,
            rename=True
        ),
    )
    process_item.var_formulas = var_formulas
    process_item.bounds = bounds
    try:
        # initialize template function
        process_item.template = None
        if template_file:
            process_item.template = compile_template(template_file)

    except Exception as e:
        if cli_arguments.debug:
            traceback.print_exc()
        elif isinstance(e, TemplateError):
            traceback.print_exc(limit=0)
        exit_program("Error in template file: " + str(e))


def init_pool_process(
        par_names, var_names, data_names, file_column_names,
        var_formulas, bounds, template_file, config_dir, helper_modules
):
    """Sub-process version of init_processing.

    Also switches off SIGINT. which is handled by mother process.
    """
    # Ignore SIGINT (KeyboardInterrupt) and warnings if started as pool worker
    # process. Only the mother process handles these and terminates children

    signal.signal(signal.SIGINT, signal.SIG_IGN)
    warnings.filterwarnings("ignore")
    # do not show code in warnings
    warnings.showwarning = lambda message, category, filename, lineno, \
        file=sys.stderr, line=None, original=warnings.showwarning: \
        original(message, category, filename, lineno, file, "")

    init_processing(
        par_names, var_names, data_names, file_column_names,
        var_formulas, bounds, template_file, config_dir, helper_modules
    )


class TemplateError(SyntaxError):

    """Exception for syntax errors in templates."""

    pass


def compile_template(template_file):
    """Return compiled code object for given template.

    Code object evaluates to template string with placeholders substituted
        with values in local and global scope

    Allowed placeholders:
        $$ -> evaluated to "$"
        $identifier -> evaluated to repr(identifier)
        ${extidentifier} -> evaluated to repr(extidentifier)

        where identifier is a string starting with "_" or a letter and
            continuing with letters, numbers, "." or "_"
            and extidentifier is the same but can also contain "[]"
            characters in the continuing part

    Note: this function uses process_item.with_names, so it can only be
        called after or at the end of init_processing
    """
    if os.path.isfile(template_file):
        with open(template_file) as f:
            template = f.read()
    else:
        template = template_file
        template_file = "<template string>"

    (
        pars_with_names, vars_with_names, par_and_vars_with_names
    ) = process_item.with_names
    pars_example = pars_with_names(*([0] * len(pars_with_names._fields)))
    vars_example = vars_with_names(*([0] * len(vars_with_names._fields)))
    values_example = par_and_vars_with_names(
        *([0] * len(par_and_vars_with_names._fields))
    )

    pattern = r"""\$(?:
        (\$)|                       # escape sequence
        ([_a-z][_a-z0-9.]*)|        # simple identifier
        {([_a-z][_a-z0-9.\[\]]*)}|  # extended identifier
        ()                          # error
    )"""
    pos = 0
    parts = []
    error = None
    for match in re.finditer(pattern, template, re.IGNORECASE + re.VERBOSE):
        # group 1: $$
        # group 2: $identifier
        # group 3: ${extidentifier}
        # group 4: error
        parts.append(repr(template[pos:match.start()]))
        key = match.group(2) or match.group(3)
        if match.group(1) is not None:
            parts.append(repr("$"))
        elif key is not None:
            parts.append("repr(" + key + ")")
            try:
                bytecode = compile(key, template_file, "eval")
                eval(
                    bytecode,
                    {
                        "pars": pars_example,
                        "vars": vars_example,
                        "values": values_example
                    },
                    values_example._asdict()
                )
            except SyntaxError as e:
                # add 1 to offset to accommodate for "{" char in case of
                #   braced placeholder
                if match.group(3) is not None:
                    e.offset += 1
                error = '%s: %s' % (e.msg, e.text), match.start() + e.offset
                break
            except (AttributeError, NameError, IndexError) as e:
                error = str(e), match.start()
                break
        elif match.group(4) is not None:
            error = "invalid placeholder", match.start()
            break
        else:
            error = 'unrecognized group in pattern', match.start()
            break
        pos = match.end()

    if error is not None:
        error_str, error_loc = error
        lineno = 1
        col = 1
        line = ""
        for c in template[:error_loc]:
            if c == "\n":
                lineno += 1
                col = 1
                line = ""
            else:
                col += 1
                line += c
        for c in template[error_loc:]:
            if c == "\n":
                break
            line += c

        raise TemplateError(error_str, (template_file, lineno, col, line))

    parts.append(repr(template[pos:]))
    parts = [p for p in parts if p != repr("")]

    # building the string using ''.join is about 30% faster than using the
    #   string addition operator
    return compile(
        "''.join([" + ', \n'.join(parts) + "])\n",
        "code:" + template_file,
        "eval"
    )


def color(c, text):
    """Return text with VT100 colors."""
    if c == "red":
        return "\x1B[31m" + text + "\x1B[0m"
    if c == "green":
        return "\x1B[32m" + text + "\x1B[0m"
    if c == "yellow":
        return "\x1B[33m" + text + "\x1B[0m"
    return text


def exit_program(message=None, successful=False):
    """Exit program with message to terminal/log file.

    Also cleans up multiprocessing pool.
    """
    if message is not None:
        warn(message, seriousness=-1 if successful else 1)

    if hasattr(process_batch, 'pool'):
        process_batch.pool.terminate()
        process_batch.pool.join()

    sys.exit(0 if successful else 1)


def warn(*message, **kwargs):
    """Print warning message to stderr (in color if supported).

    Keyword arguments:
        seriousness     can be -1(informational), 0(warning), 1(error)
                            defaults to 0
        log             bool, controls writing message also to file
                            defaults to True
    """
    message = ' '.join(map(str, message))
    seriousness = kwargs.pop("seriousness", 0)
    log = kwargs.pop("log", True)
    if kwargs:
        raise TypeError(
            "warn does not take these arguments: " +
            ', '.join(kwargs)
        )

    if not message.endswith("\n"):
        message += "\n"

    # log before adding color
    if log and hasattr(warn, "log_file"):
        with open(warn.log_file, "a") as f:
            f.write(message)

    # warning messages are in yellow (if possible) or red for serious stuff
    # negative seriousness indicates success so green
    if sys.stderr.isatty():
        message = color(
            {-1: "green", 0: "yellow", 1: "red"}[seriousness],
            message
        )

    sys.stderr.write(message)
    sys.stderr.flush()


def error(msg):
    """Shortcut to warn(*, seriousness=1) for error messages."""
    warn(msg, seriousness=1)


def log_msg(message):
    """Write message to log file."""
    if hasattr(warn, "log_file"):
        with open(warn.log_file, "a") as f:
            f.write(message)
            f.write("\n")


# ----------------------------------- #
# | [E] Configuration File Handling | #
# ----------------------------------- #

# -------------------------------- #
# |  [E.1] Processing Directives | #
# -------------------------------- #

# set some reasonable defaults
config = ConfigParser.SafeConfigParser()
config.optionxform = str  # disable automatic lower case'ing
config.add_section("parameter_space")
config.set("parameter_space", "bound_count", "0")
config.add_section("algorithm")
config.set("algorithm", "projection_count", "0")
config.set("algorithm", "symmetry_count", "0")

no_config = False
if not cli_arguments.input_file:
    no_config = True
else:
    chosen_scan_setup = cli_arguments.input_file[-1].name
    for in_file in cli_arguments.input_file:
        file_name = in_file.name
        # we read it through the preprocessor, so close this file instance
        in_file.close()
        try:
            with PreprocessedConfigFile(file_name) as f:
                config.readfp(f, filename=file_name)
        except IOError as e:
            exit_program("Error: " + str(e))
        except ConfigParser.ParsingError as e:
            exit_program("Error: " + str(e))

# no config file given or config file not readable
#   -> let user choose one (if there is one in this folder)
if no_config:
    if not sys.stdin.isatty():
        exit_program("No scan definition to process")
    scan_setups = glob("*.scan")
    if not scan_setups:
        exit_program("No scan definition files found")

    scan_setups.sort()
    scan_setup_names = []

    # NOTE the list is 1-based, but Python arrays are 0-based
    for (i, m) in enumerate(scan_setups):
        scan_setup_names.append(m)
        print "[", (i + 1), "]", scan_setup_names[-1]

    chosen_scan_setup = None
    while chosen_scan_setup is None:
        try:
            # default value (if user just presses Enter) is the first found
            #   config file (in the glob list)
            choice = raw_input("# Which model? [1] ").strip() or "1"
            if choice.isdigit():
                m = int(choice) - 1
                chosen_scan_setup = scan_setups[m]
                break
            chosen_scan_setup = choice
        except (KeyboardInterrupt, EOFError) as e:
            print  # print one last line for cosmetic purposes
            exit_program()
        except (ValueError, IndexError) as e:
            exit_program("Error: %s" % e)
        except Exception as e:
            # Ask until the value fits
            exit_program("Unexpected Error: %s" % e)

    try:
        with PreprocessedConfigFile(chosen_scan_setup) as f:
            config.readfp(f, filename=chosen_scan_setup)
    except (IOError, ConfigParser.ParsingError) as e:
        exit_program("Error: " + str(e))

config_dir = os.path.abspath(os.path.dirname(chosen_scan_setup))
print('config_dir is: ', config_dir)
base_name = os.path.basename(chosen_scan_setup)
print('base_name is: ', base_name)
draw_progress.scan_file = os.path.abspath(chosen_scan_setup)

if cli_arguments.output_dir:
    try:
        if not os.path.isdir(cli_arguments.output_dir):
            os.makedirs(
                os.path.abspath(cli_arguments.output_dir)
            )
            print('os.path.abspath(cli_arguments.output_dir: ', os.path.abspath(cli_arguments.output_dir))
    except Exception as e:
        exit_program("Error: %s" % e)

    base_name = os.path.join(cli_arguments.output_dir, base_name)
    print("base_name is: ", base_name)
warn.log_file = base_name + ".log"
warn("# Starting at", time.strftime("%d.%m.%Y %H:%M"))
atexit.register(
    lambda: warn("# Stopping at", time.strftime("%d.%m.%Y %H:%M"))
)


def config_get(type, section, name,
               default=Exception, config=config, raw=False):
    """Helper function for reading directives (including error handling)."""
    try:
        result = {
            int: config.getint,
            float: config.getfloat,
            str: (lambda s, k: config.get(s, k, raw=raw))
        }[type](section, name)
        config_get.processed_options.add((section, name))
        return result
    except (ConfigParser.NoSectionError, ConfigParser.NoOptionError,
            ConfigParser.InterpolationError) as e:
        if default is not Exception:
            return default
        exit_program("Error: %s" % e)
    except ValueError as e:
        if default is not Exception:
            print "Warning for [%s] %s (defaulting to %r): %s" % (
                section, name, default, e
            )
            return default
        exit_program("Error in [%s] %s: %s" % (section, name, e))


config_get.processed_options = set()
# intercept all ConfigParser.get calls to keep track of processed options
config.get = lambda section, option, raw=False, vars=None, getter=config.get: [
    config_get.processed_options.add((section, option)),
    getter(section, option, raw, vars)
][1]

print "# Processing", chosen_scan_setup
warn("# Results will be saved to:", os.path.abspath(base_name + ".*"))

# ------------------------------- #
# |  [E.1.1] General Directives | #
# ------------------------------- #

# config files have to contain "setup" and "parameter_space" blocks and the
#   setup block has to specify the point processor
if not config.has_section("setup"):
    exit_program("Incomplete setup given: setup section missing")

# furthermore, the version of the config file has to match the program version
scan_version = config_get(int, "setup", "version", default=1)
if not 1 <= scan_version <= __major_version__:
    exit_program(
        "Error: incompatible version, got %i, expecting between 1 and %i" % (
            scan_version, __major_version__
        )
    )
elif scan_version != __major_version__:
    print "# Running scan backwards-compatible to", scan_version

if not config.has_section("parameter_space"):
    exit_program("Error: incomplete setup given")

# default mode is "scan" (also fallback if unknown mode is given)
if cli_arguments.mode is not None:
    # cli_arguments.mode is always constrained to be a key to mode_names by
    #   the argparse module
    mode = mode_names[cli_arguments.mode]
elif config.has_option("setup", "mode"):
    mode_name = config_get(str, "setup", "mode").strip().lower()
    try:
        mode = mode_names[mode_name]
    except KeyError:
        exit_program("Error: unrecognized mode name '%s'" % mode_name)
else:
    mode = MODE_SCAN

# seed for random number generator
seed = cli_arguments.randomseed or hash(random.SystemRandom().random())
warn("# Using random seed: %r" % seed)
random.seed(seed)

concurrent_processors = config_get(
    int, "setup", "concurrent_processors",
    default=multiprocessing.cpu_count()
)
processing_concurrency.value = concurrent_processors

# parspace mode specifies whether the points to be calculated
#   are given as parameter ranges ("grid") or
#   are taken from a predetermined file containing coordinates for them
#       ("file")
if config.has_option("parameter_space", "mode"):
    parspace_mode = config_get(str, "parameter_space", "mode").lower()
    if parspace_mode not in parspace_mode_names:
        exit_program(
            "Error: unrecognized parameter space mode: %s" % parspace_mode
        )
    parspace_mode = parspace_mode_names[parspace_mode]
else:
    parspace_mode = PARSPACE_MODE_GRID

# this directive is only available for "scan" and "explore" modes
parspace_files = []
if (mode == MODE_SCAN and parspace_mode == PARSPACE_MODE_FILE) or \
   (mode == MODE_EXPLORER and config.has_option("parameter_space", "files")):
    # parspace files are taken from one config option "files" in
    #   parameter_space block and are supposed to be separated by ":"
    parspace_files = []
    for file_name in config_get(str, "parameter_space", "files").split(":"):
        if cli_arguments.output_dir and not os.path.dirname(file_name) and \
           os.path.isfile(os.path.join(cli_arguments.output_dir, file_name)):
            file_name = os.path.join(cli_arguments.output_dir, file_name)
        else:
            file_name = find_file(file_name, config_dir)

        if not os.path.isfile(file_name):
            exit_program("Error: not a file: " + file_name)

        parspace_files.append(file_name)

if mode not in [MODE_SCAN, MODE_TEST] and \
   parspace_mode != PARSPACE_MODE_GRID:
    warn(
        "# Warning: this parameter space mode is not supported in this mode: "
        +
        dict(
            zip(parspace_mode_names.values(), parspace_mode_names.keys())
        )[parspace_mode]
    )
    warn("# Reverting to default parameter space mode")
    parspace_mode = PARSPACE_MODE_GRID

# import the point_processor module and initialize it if necessary
#   this can e.g. search for binary files used for the actual calculations
# the point_processor module must have the function "main" and can have
#   the function "init"
point_processor = config_get(str, "setup", "point_processor", default=None)
if point_processor:
    point_processor = find_file(
        point_processor,
        config_dir
    )
    if not os.path.isfile(point_processor):
        exit_program("Error: not a file: " + point_processor)

    # any errors generated by this are NOT caught by the main program and are
    #   shown to the user without any filtering
    point_processor = imp.load_source("processor_module", point_processor)
    if "init" in dir(point_processor):
        point_processor.init(config_dir, config, sys.modules[__name__])
else:
    warn("# Warning: empty `point_processor' directive")

template_file = config_get(str, "setup", "template", default=None)
if template_file:
    template_file = find_file(template_file, config_dir)
    if not os.path.isfile(template_file):
        exit_program("Error: not a file: " + template_file)
else:
    warn("# Warning: empty `template' directive")

if config.has_option("setup", "helper_modules"):
    helper_modules = map(
        str.strip,
        config_get(str, "setup", "helper_modules").split(":")
    )
    helper_modules = filter(bool, helper_modules)
else:
    helper_modules = []

par_names = map(
    str.strip, config_get(str, "parameter_space", "par_names").split(",")
)
par_count = len(par_names)
if not par_count:
    exit_program("Error: no parameters defined")
try:
    check_names(par_names)
except Exception as e:
    exit_program("Error in par_names: %s" % e)

# Beware of ''.split(",") == ['']!
var_names = config_get(str, "parameter_space", "var_names", default="")
var_names = map(str.strip, var_names.split(",")) if var_names else []
var_count = len(var_names)
try:
    check_names(var_names)
except Exception as e:
    exit_program("Error in var_names: %s" % e)

data_names = config_get(str, "parameter_space", "data_names", default="")
data_names = map(str.strip, data_names.split(",")) if data_names else []
try:
    check_names(data_names)
except Exception as e:
    exit_program("Error in data_names: %s" % e)

if data_names and not point_processor:
    exit_program(
        "Error: empty `point_processor' directive,"
        " but nonempty `data_names'"
    )

bound_count = config_get(int, "parameter_space", "bound_count")

if mode in [MODE_MCMC, MODE_EXPLORER, MODE_OPTIMIZE] and \
   not config.has_option("algorithm", "likelihood"):
    # exit program if likelihood missing in mode that needs it
    exit_program("Error: no likelihood specified")

if mode in [MODE_MCMC, MODE_EXPLORER, MODE_OPTIMIZE] and \
   config.has_option("algorithm", "out_columns"):
    # out_columns is only used in test and scan mode
    #   in worker mode, this is ignored entirely
    warn("# Warning: out_columns only supported in scan mode")

# access columns of input files directly
if parspace_mode == PARSPACE_MODE_FILE and \
   config.has_option("parameter_space", "file_columns"):
    file_column_names = map(
        str.strip,
        config_get(str, "parameter_space", "file_columns").split(",")
    )
else:
    file_column_names = []

par_definitions = [
    config_get(str, "parameter_space", "par_%s" % n)
    for n in par_names
]
raw_var_formulas = [
    config_get(str, "parameter_space", "var_%s" % n, raw=True)
    for n in var_names
]
var_formulas = [
    config_get(str, "parameter_space", "var_%s" % n)
    for n in var_names
]
bounds = [
    config_get(str, "parameter_space", "bound_%i" % i)
    for i in range(bound_count)
]

# | End of general configuration
# |----------------------------------------------------------------------------
# | Next up mode specific configuration (after creating pool)

# this initializes formula_eval and process_item
try:
    init_processing(
        par_names, var_names, data_names, file_column_names,
        var_formulas, bounds, template_file, config_dir, helper_modules
    )
except Exception as e:
    if cli_arguments.debug:
        traceback.print_exc()
    exit_program("Error during initialization: " + str(e))

# create the pool now, so the children have their own copies of all parameter
#   space info but not of the scan specific data
# could be especially costly for explorer mode
# also test mode does not use the (local) pool at all
if processing_concurrency() and mode != MODE_TEST:
    process_batch.pool = multiprocessing.Pool(
        processes=processing_concurrency(),
        initializer=init_pool_process,
        initargs=(
            par_names, var_names, data_names, file_column_names,
            var_formulas, bounds, template_file, config_dir, helper_modules
        )
    )
else:
    # also gracefully handle situation of no local calculation allowed
    #   i.e. have pool object that does nothing without errors
    class inert(object):

        """Object that does nothing and raises no errors doing it."""

        def __init__(self):
            """Do nothing."""
            pass

        def __getattribute__(self, name):
            """Return only None as attribute."""
            return lambda: None

    process_batch.pool = inert()

# -----------------------------------------------------------------------------

# ------------------------------ #
# |  [E.1.2] Shared Directives | #
# ------------------------------ #

# check validity of parameter ranges
if parspace_mode != PARSPACE_MODE_FILE:
    par_ranges = []
    effective_par_count = 0
    for i in range(par_count):
        try:
            par_ranges.append(ParameterRange(par_definitions[i]))
            if not (par_ranges[i].is_finite and len(par_ranges[i]) == 1):
                effective_par_count += 1
        except ValueError as e:
            exit_program(
                "Parameter range %r invalid: %s" % (par_definitions[i], e)
            )
else:
    # in PARSPACE_MODE_FILE, the setup comes later
    par_ranges = []
    effective_par_count = par_count
    for i in xrange(par_count):
        try:
            compile(
                "(" + par_definitions[i] + ")",
                abbreviate(repr(par_definitions[i])),
                "eval"
            )
        except Exception as e:
            exit_program("Error in formula for parameter %i: %s" % (i, e))

# unit_length has different interpretations depending on mode
# scan: number of points that are calculated between saving
# mcmc: length of the chains
# explore: number of points that are simultaneously explored
# optimize: size of population
# test and worker do not use it
default_unit_length = lambda: {
    MODE_SCAN: 100 * concurrent_processors,
    MODE_OPTIMIZE: max(10 * effective_par_count, 1),
    MODE_EXPLORER: max(10 * effective_par_count, 1),
    MODE_TEST: 0,
}
unit_length = config_get(
    int, "setup", "unit_length",
    default=None if mode in default_unit_length() else Exception
)
unit_length_default_used = False
if unit_length is None:
    unit_length_default_used = True
    unit_length = default_unit_length()[mode]

# check validity of var and bound formulas (only syntax)
has_errors = False
for i, formula in enumerate(var_formulas):
    try:
        compile("(" + formula + ")", abbreviate(repr(formula)), "eval")
    except Exception as e:
        error("Error in formula for %r: %s" % (var_names[i], e))
        has_errors = True

for i, formula in enumerate(bounds):
    try:
        compile("(" + formula + ")", abbreviate(repr(formula)), "eval")
    except Exception as e:
        error("Error in formula for bound %i: %s" % (i, e))
        has_errors = True

if has_errors:
    exit_program()

if mode in [MODE_MCMC, MODE_EXPLORER, MODE_OPTIMIZE] or \
   (config.has_option("algorithm", "likelihood") and mode == MODE_TEST):
    # only define likelihood function if needed
    likelihood_code = config_get(str, "algorithm", "likelihood")
    try:
        compile(
            "(" + likelihood_code + ")",
            abbreviate(repr(likelihood_code)),
            "eval"
        )
    except Exception as e:
        exit_program("Error in formula for likelihood: " + str(e))

    def likelihood(data, code=likelihood_code):
        """Return likelihood of given parameter point data triple."""
        return formula_eval(
            code,
            pars=data[0],
            vars=data[1],
            data=data[2]
        )

else:
    likelihood_code = None

# only used as reference for test mode; for explorer mode, see below
min_likelihood = config_get(
    float,
    "algorithm", "min_likelihood",
    default=None if mode == MODE_TEST else Exception
) if mode in [MODE_TEST, MODE_EXPLORER] else None

# -----------------------------------------------------------------------------
# | projections and symmetries are only used by explorer mode, but also shown |
# |  in test mode for checking                                                |
# -----------------------------------------------------------------------------

# due to setting a default value of 0 previously, the default value given
#   here only applies in cases of ValueError and
#   ConfigParser.InterpolationError etc.
projection_count = config_get(
    int, "algorithm", "projection_count",
    default=0 if mode == MODE_TEST else Exception
)
projections = []
has_errors = False
for i in range(projection_count):
    # single coordinate form takes precedence over combined one, but only
    #   if both are given
    if not config.has_option("algorithm", "projection_%i_x" % i) or \
       not config.has_option("algorithm", "projection_%i_y" % i):

        if not config.has_option("algorithm", "projection_%i" % i):
            exit_program(
                "Error: incomplete projection specification for "
                "projection %i" % i
            )
        parts = split_by_comma(
            config_get(str, "algorithm", "projection_%i" % i)
        )
        if len(parts) != 2:
            exit_program(
                "Error: malformed projection specification for projection "
                "%i: expected 2 coordinates (%i given)" % (i, len(parts))
            )

        x_formula, y_formula = map(str.strip, parts)
    else:
        x_formula = config_get(str, "algorithm", "projection_%i_x" % i).strip()
        y_formula = config_get(str, "algorithm", "projection_%i_y" % i).strip()

    if config.has_option("algorithm", "projection_%i_z" % i):
        z_formula = config_get(str, "algorithm", "projection_%i_z" % i).strip()
    else:
        z_formula = "likelihood"

    if config.has_option("algorithm", "projection_%i_filter" % i):
        filter_formula = config_get(
            str, "algorithm", "projection_%i_filter" % i
        ).strip()
    else:
        filter_formula = ""

    projections.append((x_formula, y_formula, z_formula, filter_formula))

    # check each coordinate for syntax errors
    named_formulas = zip(
        ("x coordinate", "y coordinate", "z coordinate", "filter"),
        projections[-1]
    )
    for name, formula in named_formulas:
        try:
            if (name != "filter" or formula) and formula != "likelihood":
                compile(
                    "(" + formula + ")",
                    abbreviate(repr(formula)),
                    "eval"
                )
        except Exception as e:
            error("Error in formula for %s: %s" % (name, e))
            has_errors = True

# due to setting a default value of 0 previously, the default value given
#   here only applies in cases of ValueError and
#   ConfigParser.InterpolationError etc.
symmetry_count = config_get(
    int, "algorithm", "symmetry_count",
    default=0 if mode == MODE_TEST else Exception
)
symmetry_trafos = []
for i in range(symmetry_count):
    try:
        code = config_get(str, "algorithm", "symmetry_%i" % i)
        compile(
            "{ %s }" % code,
            abbreviate(repr(code)),
            "eval"
        )
        symmetry_trafos.append("{ %s }" % code)
    except Exception as e:
        error("Error in formula for symmetry %i: %s" % (i, str(e)))
        has_errors = True

if has_errors:
    exit_program()

# -------------------------------- #
# |  [E.1.3] Specific Directives | #
# -------------------------------- #

out_columns = None
if mode == MODE_SCAN or mode == MODE_TEST:
    # "out_columns" can be used as selector (and transformation) of which
    #   data should be written as result for each point
    out_columns = config_get(str, "algorithm", "out_columns", default=None)
    if out_columns is not None:
        out_columns = map(str.strip, split_by_comma(out_columns))
        for i, col in enumerate(out_columns):
            try:
                compile(
                    "(" + col + ")",
                    abbreviate(repr(col)),
                    "eval"
                )
            except Exception as e:
                exit_program(
                    "Error in formula for output column %i: %s" % (i, e)
                )

        # define post_processor, which selects/transforms output columns
        #   according to config and returns a single string array
        def post_processor(pars, vars, data, file=None, columns=out_columns):
            """Return point values as requested in out_columns directive."""
            return [
                out_repr(
                    formula_eval(
                        col, pars=pars, vars=vars, data=data, file=file
                    )
                )
                for col in columns
            ]
    else:
        def post_processor(pars, vars, data, file=None):
            """Return point values as per default."""
#           Warning: Changed here
#           return map(out_repr, pars + vars + data)
            return map(out_repr, data)

    if mode == MODE_SCAN and parspace_mode == PARSPACE_MODE_SCATTER:
        scatter_point_count = config_get(int, "parameter_space", "point_count")

elif mode == MODE_MCMC:
    heating_function_code = config_get(
        str, "algorithm", "heating_function", default=None
    )
    if heating_function_code:
        try:
            compile(
                "(" + heating_function_code + ")",
                abbreviate(repr(heating_function_code)),
                "eval"
            )
        except Exception as e:
            exit_program("Error in formula for heating function: " + str(e))

elif mode == MODE_OPTIMIZE:
    differential_weight = config_get(
        float, "algorithm", "differential_weight", default=0.6
    )
    crossover_probability = config_get(
        float, "algorithm", "crossover_probability", default=0.5
    )
    epsilon = abs(
        config_get(float, "algorithm", "waiting_threshold", default=0)
    )
    # set default requested precision to half machine precision
    epsilon_relative = abs(
        config_get(
            float, "algorithm", "waiting_threshold_relative",
            default=10 ** round(math.log10(sys.float_info.epsilon) / 2)
        )
    )
    max_waiting_time = config_get(
        int, "algorithm", "waiting_time",
        default=10 * effective_par_count
    )
    if unit_length < 4:
        exit_program(
            ("Error: population size '%i' too small " % unit_length) +
            "(has to be larger than 4)"
        )

elif mode == MODE_EXPLORER:
    # number of points which are used as exploration basis each iteration
    # restriction does not apply to interpolated and extrapolated points
    unit_length  # (already processed)

    # minimal likelihood considered
    # points with less likelihood are ignored and treated as invalid (but
    #   will show up in data file and not excluded data file)
    min_likelihood  # (already processed)

    # breaks up points into categories depending on likelihood
    #   points in lower likelihood bin will only be considered for
    #   exploration if bins with higher likelihood have been depleted
    if config.has_option("algorithm", "likelihood_steps"):
        # sort the likelihood bins starting with the highest in descending
        #   order
        likelihood_steps = sorted(
            set(
                map(
                    float,
                    config_get(str, "algorithm", "likelihood_steps").split(",")
                ) + [min_likelihood]
            ),
            reverse=True
        )
        likelihood_steps = tuple(likelihood_steps)
    else:
        likelihood_steps = (min_likelihood, )

    if config.has_option("algorithm", "disabled_states"):
        available_states = set([1, 2, 3]) - set(map(
            int,
            filter(
                bool,
                config_get(
                    str, "algorithm", "disabled_states"
                ).split(",")
            )
        ))
        available_states = sorted(available_states)
        if not available_states:
            exit_program("Error: all algorithm states disabled")
        if not projections and 3 not in available_states:
            exit_program(
                "Error: no projections defined and algorithm state 3 disabled"
            )
        if not projections:
            # state one and two only work with projections
            available_states = sorted(set(available_states) - set([1, 2]))
    else:
        available_states = [1, 2, 3]
        if not projections:
            available_states = sorted(set(available_states) - set([1, 2]))

    # purge points from the cache if their likelihood falls into a bin that
    #   has long been depleted
    min_purge_age = config_get(int, "algorithm", "min_purge_age", default=2)
    if min_purge_age <= 0:
        # never purge the source_bin or worse
        exit_program(
            "Error: minimal purge age must be larger than 0: %i" %
            min_purge_age
        )

    # only load points that pass the specified loading filter
    if config.has_option("algorithm", "loading_filter"):
        loading_filter_code = config_get(str, "algorithm", "loading_filter")
        try:
            compile(
                "(" + loading_filter_code + ")",
                abbreviate(repr(loading_filter_code)),
                "eval"
            )
        except Exception as e:
            exit_program("Error in formula for loading filter: " + str(e))

        loading_filter = lambda data, code=loading_filter_code: formula_eval(
            code,
            pars=data[0],
            vars=data[1],
            data=data[2]
        )
    else:
        loading_filter = lambda data: True

    # if extrapolation is directly specified for only some projections,
    #   restrict to those, otherwise extrapolate on all projections
    if config.has_option("algorithm", "extrapolated_projections"):
        try:
            # filter to handle value '' correctly (it gets split to ['']
            extrapolated_projections = map(
                int,
                filter(
                    bool,
                    config_get(
                        str, "algorithm", "extrapolated_projections"
                    ).split(",")
                )
            )
        except Exception as e:  # most certainly ValueError from int()
            exit_program("Error: " + str(e))
    else:
        extrapolated_projections = range(projection_count)

    # use suspected symmetries in the data to fill up parameter space evenly
    #   i.e. for every point also its symmetry-partner points are queued for
    #   calculation only works with parameters accessed via dictionary "pars"
    # format: list of "input parameter name: output parameter value", where
    #   the output value is a function of the input values stored in the
    #   namedtuple "pars"
    symmetry_count   # (already processed)
    symmetry_trafos  # (already processed)

if config.has_option("setup", "workers") and mode != MODE_MCMC:
    # configuration for manager mode
    authkey = config_get(str, "setup", "authkey")

    # disable workers if in worker mode or no workers are specified
    worker_addresses = None
    if not worker_addresses:
        worker_addresses = None
else:
    authkey = None
    worker_addresses = None

# Show the user what other directives were present but ignored

# this is set by default but not always used, so "take" it manually
#   (but raise no exception)
config_get(int, "algorithm", "symmetry_count", default=0)

ignored = []
for sec in config.sections():
    for k, v in config.items(sec):
        if (sec, k) not in config_get.processed_options:
            ignored.append((sec, k))
if ignored:
    warn(
        "# Warning: the following directives have been ignored "
        "(not counting text interpolation):"
    )
    for s_and_k in sorted(ignored):
        warn("\t[%s] %s" % s_and_k)

# ------------------------- #
# |  [E.2] Review by user | #
# ------------------------- #

# no access to config.* and config_get after this point!
config = None
config_get = None

print "# Please review the scan definition " \
    "and check that all formulas are safe to run:"

indent = "    "

warn("# Mode:", dict(zip(mode_names.values(), mode_names.keys()))[mode])
warn("# Concurrent processors (locally):", processing_concurrency())
if mode not in [MODE_TEST]:
    warn("# Unit length:", unit_length)

if parspace_files:
    print "# Data read from files:"
    for file_name in parspace_files:
        print file_name
    if file_column_names:
        print "# in format:", ", ".join(file_column_names)

# parspace_mode == PARSPACE_MODE_FILE is only available for scan
#   everything else is PARSPACE_MODE_GRID
print "# Parameters:"
if par_count:
    for i in range(par_count):
        print "pars[%i] =" % i, par_names[i], "=",  # note dangling ','
        if parspace_mode == PARSPACE_MODE_GRID:
            print par_ranges[i]
            if cli_arguments.debug:
                print textwrap.fill(
                    repr(par_ranges[i].values),
                    initial_indent=indent,
                    subsequent_indent=indent,
                    width=terminal_size()[0]
                )
                print "%s%i values" % (
                    indent, len(par_ranges[i].values)
                )
        else:
            print par_definitions[i]

else:
    print "-- none --"

print "# Variables (interpolation not shown):"
if var_count:
    for i in range(var_count):
        # indent new lines in formulas
        print "vars[%i] =" % i, var_names[i], "=",  # note dangling ','
        print raw_var_formulas[i].replace("\n", "\n" + indent)
else:
    print "-- none --"

if likelihood_code is not None:
    print "# Likelihood function:"
    print "likelihood =", likelihood_code.replace("\n", "\n" + indent)
    if min_likelihood is not None and mode != MODE_EXPLORER:
        print "min_likelihood =", min_likelihood

    print "# Output columns are then (separated by \"\\t\", " \
        "shown as \"index : name\", possibly incomplete):"

out_columns_for_review = par_names + var_names + data_names
if mode == MODE_MCMC:
    out_columns_for_review += ["likelihood", "staycount"]
elif mode in [MODE_OPTIMIZE, MODE_EXPLORER]:
    out_columns_for_review += ["likelihood"]
elif mode in [MODE_SCAN, MODE_TEST]:
    if out_columns is not None:
        out_columns_for_review = out_columns
    else:
        out_columns_for_review = par_names + var_names + data_names
else:
    # really raise an uncaught exception since this is a coding error
    raise ValueError("Unknown mode")

print_table(
    map(str, range(1, len(out_columns_for_review) + 1)),
    out_columns_for_review
)

if out_columns is not None and cli_arguments.debug:
    print "# Raw output columns: ", out_columns

print "# Bounds:"
if bound_count:
    for i, bound in enumerate(bounds):
        # indent new lines in formulas
        print "bound_%i =" % i, bound.replace("\n", "\n" + indent)
else:
    print "-- none --"

if mode in [MODE_TEST, MODE_EXPLORER]:
    print "# Projections:"
    if projection_count:
        for proj in projections:
            x_formula, y_formula, z_formula, filter_formula = proj
            if filter_formula:
                print "(%s, %s, %s) if %s" % (
                    x_formula, y_formula, z_formula, filter_formula
                )
            else:
                print "%s, %s, %s" % (x_formula, y_formula, z_formula)
    else:
        print "-- none --"

if mode == MODE_OPTIMIZE:
    print "# Algorithm specific options:"
    print "#    Differential weight:", differential_weight
    print "#    Crossover probability:", crossover_probability
    print "#    Likelihood convergence threshold:", epsilon
    print "#    Likelihood convergence threshold (relative):", epsilon_relative
    print "#    Convergence waiting time:", max_waiting_time
elif mode == MODE_EXPLORER:
    print "# Algorithm specific options:"
    print "#    Minimal likelihood for further consideration:", min_likelihood
    if len(likelihood_steps) > 1:
        print "#    Boundary likelihood bins:", ', '.join(
            [
                "[%r, %r)" % (
                    likelihood_steps[i],
                    likelihood_steps[i - 1] if i - 1 >= 0 else "inf"
                )
                for i in range(len(likelihood_steps) - 1, -1, -1)
            ]
        )
    if projection_count:
        print "#    Extrapolated projections:", extrapolated_projections
elif mode == MODE_MCMC:
    if heating_function_code:
        print "# Algorithm specific options:"
        print "#    Heating function:", heating_function_code

# ------------------------- #
# | [F] Manager Mode Setup | #
# ------------------------- #

# setup the appropriate calculate_items function: dispatches jobs to all
#   available workers/process and returns iterator (or iterable object) over
#   the results
if worker_addresses is None:
    # Run calculation only locally.
    calculate_items = process_batch
else:
    print "# Starting in manager mode"
    # register functions for (remote) use
    BaseManager.register("process_batch")
    BaseManager.register("process_item")
    BaseManager.register("processing_concurrency")
    RemoteMachine = namedtuple("RemoteMachine", ["addr", "concurrency"])

    # the local computer is also used for calculations if specified
    workers = [RemoteMachine(None, processing_concurrency())]
    # connect to all workers given in the comma delimited list specified by
    #  the "workers" directive in the setup block
    # every worker has to have the form "hostname:port" with an integer port
    # get the worker's processing power and disconnect again
    #   (workers are reconnected on every batch of items)
    for worker_desc in worker_addresses.split(","):
        try:
            worker = None
            host, colon, port = worker_desc.strip().partition(':')
            port = int(port) if port else DEFAULT_PORT
            worker_addr = (host, port)
            worker = BaseManager(address=worker_addr, authkey=authkey)
            worker.connect()
            worker_concurrency = worker.processing_concurrency()._getvalue()
        except (IOError, multiprocessing.connection.AuthenticationError) as e:
            warn("Could not connect to '%s': %s" % (worker_desc, e))
        except Exception as e:
            traceback.print_exc()
            warn("Could not connect to '%s': %s" % (worker_desc, e))
        else:
            # only keep if there was no error and already added
            machine = RemoteMachine(worker_addr, worker_concurrency)
            if machine not in workers:
                workers.append(machine)
        finally:
            # close all connections in any case
            del worker

    # from here on concurrent_processors also takes into account processes
    #   that live in pools of worker processes
    concurrent_processors = sum(worker.concurrency for worker in workers)

    # weigh workers by their relative available computing power
    # this means that if concurrent_processors is 0 in the manager scan file,
    #   the local manager computer will not be attributed any work
    warn(
        "# Calculation identified with authkey '%s' distributed over:" %
        authkey
    )
    warn("#    local (concurrency=%s)" % workers[0].concurrency)
    for worker in workers[1:]:
        warn(
            "#    %s (weight=%s)" % (
                ":".join(map(str, worker.addr)),
                worker.concurrency
            )
        )
    warn("# Total number of processors: ", concurrent_processors)
    if unit_length_default_used and concurrent_processors > \
       processing_concurrency():
        unit_length = default_unit_length()[mode]
        warn("# Re-set unit length to:", unit_length)
    workers = tuple(workers)

    def calculate_items(items, workers=workers):
        """Distribute calculation over more than one scanner process pool."""
        if not items:
            raise StopIteration()

        concurrencies = [worker.concurrency for worker in workers]
        full_concurrency = sum(concurrencies)
        pool_iters = []
        batches = [[] for s in workers]
        distribution = []

        # we split up the requested items into chunks of size full_concurrency
        #   and cycle through all workers for that chunk with roundrobin, where
        #   each worker is used as many times as its concurrency specifies
        # example: concurrencies = [3, 2, 4]
        #   -> chunk_order = [0, 1, 2, 0, 1, 2, 0, 2, 2]
        # this is supposed to minimize waiting times due to inhomogeneous
        #   worker concurrency
        chunk_order = list(
            roundrobin(*[
                [i] * s
                for i, s in enumerate(concurrencies)
            ])
        )
        assert len(chunk_order) == sum(concurrencies)
        for i, item in enumerate(items):
            pool_idx = chunk_order[i % full_concurrency]
            batches[pool_idx].append(item)
            distribution.append(pool_idx)
        # distribution is now interleaved list of pool indices where the index
        #   specifies which is handling the item at the corresponding position

        # this means that yielding the processed items in the same order as
        #   given, this should come evenly distributed from all workers
        #   (`evenly' also taking into account relative speed differences)
        # NOTE this assumes that process_batch also maintains the order of
        #   items otherwise the algorithm will still work but not maintain
        #   order

        # generate list of iterators, each over its attributed batch of items
        pool_iters = []
        connections = []
        for worker, batch in zip(workers, batches):
            if not batch:
                pool_iters.append(None)
                continue

            if worker.addr is None:
                pool_iters.append(process_batch(batch))
            else:
                connection = BaseManager(address=worker.addr, authkey=authkey)
                connection.connect()
                pool_iters.append(connection.process_batch(batch))
                connections.append(connection)

        for idx in distribution:
            try:
                item = pool_iters[idx].next()
                yield item
            except StopIteration:
                raise RuntimeError("Internal error: premature batch end")
            except Exception as e:
                if cli_arguments.debug:
                    traceback.print_exc()
                exit_program("Error during batch calculation: " + str(e))

        # try to disconnect by deleting all connection objects
        # this may not work if something still holds references to those
        #   objects
        del connections[:]

wait_for_user()

# --------------------- #
# | [G] Scanning Code | #
# --------------------- #

if mode == MODE_TEST:
    # -------------------- #
    # |  [G.1] Test Mode | #
    # -------------------- #
    warnings.filterwarnings("always")  # don't swallow repeated warnings

    # this mode can be used to test the point_processor and workers and
    #   calculate single points on each

    # previous input point are saved between test sessions
    try:
        import readline
    except ImportError:
        try:
            import pyreadline as readline
        except ImportError:
            pass

    try:
        readline.clear_history()
        readline.read_history_file("%s.testhistory" % base_name)
    except Exception:  # can fail due to no readline or due to no history file
        pass

    import cProfile
    import pstats
    profiler = cProfile.Profile()

    interactive_mode = not bool(cli_arguments.pars)
    if interactive_mode and not sys.stdin.isatty():
        exit_program()

    test_record_separator = lambda: "%s\n%s" % (
        horizontal_line(),
        horizontal_line()
    )

    def unsaved_raw_input(prompt):
        """Return user input but don't save it in readline history."""
        answer = raw_input(prompt)
        try:
            last_saved_item = readline.get_history_item(
                readline.get_current_history_length()
            )
            if answer == last_saved_item:
                readline.remove_history_item(
                    readline.get_current_history_length() - 1
                )
        finally:
            return answer

    print test_record_separator()

    try:
        # wait for input until an exception occurs or the user interrupts
        #   with ctrl+c
        while True:
            try:
                # if we were given parameters on the command line, calculate
                #   those first
                # if we then have an active terminal, ask for parameters
                # otherwise exit
                point = None
                while cli_arguments.pars:
                    # argparse already converts to float
                    point = cli_arguments.pars.pop(0)
                    if len(point) != par_count:
                        print "Error in 'pars' argument: ",
                        print "%d values given, expected %d" % (
                            len(point), par_count
                        )
                        continue
                    break

                if point is None and not interactive_mode:
                    exit_program()

                used_random = False
                while not point:
                    try:
                        used_random = False
                        point_input = raw_input(
                            "# Enter point (format: %s; %s):\n> " % (
                                ", ".join(par_names),
                                "as numbers or 'random'"
                                if par_ranges else "as numbers"
                            )
                        )
                        if point_input.strip().lower() == "random" \
                           and par_ranges:
                            point = make_random_point(par_ranges)
                            used_random = True
                        else:
                            if "random" in point_input and par_ranges:
                                used_random = True
                            point = point_input.split(",")
                            if len(point) != par_count:
                                raise ValueError(
                                    "%d values given, expected %d" % (
                                        len(point), par_count
                                    )
                                )
                            point = [
                                make_random_point([par_ranges[i]])[0]
                                if x.strip().lower() == "random" and par_ranges
                                else float(x)
                                for i, x in enumerate(point)
                            ]
                    except ValueError as e:
                        point = None
                        print "Error:", e

                if used_random:
                    # go up one line and clear it (using VT100 codes)
                    sys.stdout.write("\x1B[1A\x1B[K")
                    print ">", ", ".join(map(repr, point))

                if worker_addresses is not None and sys.stdin.isatty():
                    prompt = "# Which worker? [%s] " % \
                        ",".join(map(str, range(len(workers))))

                    try:
                        worker_id = int(unsaved_raw_input(prompt))
                    except Exception:
                        worker_id = 0

                    worker = workers[worker_id][0]
                else:
                    worker = None

                pars_out_of_range = []
                if par_ranges:  # can be empty if launched with file mode
                    pars_out_of_range = [
                        name
                        for i, name in enumerate(par_names)
                        if point[i] not in par_ranges[i]
                    ]
                if pars_out_of_range:
                    warn(
                        "# Parameters with values not inside range or not one"
                        " of the defined values:",
                        ", ".join(pars_out_of_range),
                        log=False
                    )
                    if sys.stdin.isatty():
                        answer = unsaved_raw_input(
                            "# Do you want to round to "
                            "the nearest point? [Y/n] "
                        )
                        if answer != "n":
                            for i, x in enumerate(point):
                                point[i] = par_ranges[i].pull_inside(point[i])

                print "# Parameters:"
                print ", ".join(
                    "%s = %r" % (k, v) for (k, v) in zip(par_names, point)
                )

                if sys.stdin.isatty():
                    try:
                        readline.write_history_file(
                            "%s.testhistory" % base_name
                        )
                    except Exception:
                        pass

                starttime = time.time()
                if worker is None:
                    if cli_arguments.profiling:
                        profiler.enable()

                    data = process_item(
                        point, test_mode=True,
                        print_error_trace=cli_arguments.debug
                    )

                    # disable unconditionally so the check is not included, too
                    profiler.disable()
                    if cli_arguments.profiling:
                        ps = pstats.Stats(profiler).sort_stats(1)
                        # sort_stats(1) -> sort by "tottime"
                        ps.print_stats()
                        profiler = cProfile.Profile()
                        # recreate, otherwise each test iteration will be
                        #   appended instead of replaced
                else:
                    # remote profiling may not be such a good idea
                    worker_connection = BaseManager(
                        address=worker, authkey=authkey
                    )
                    worker_connection.connect()
                    data = worker_connection.process_batch(
                        [(point, None)]
                    ).next()[0]

                print "# Calculation done after: %f s" % (
                    time.time() - starttime
                )

                # if data is an Exception, we got an error
                if isinstance(data, Exception):
                    # don't save bad test data and just print the error message
                    data = data.args
                    print "error=\"%s\"" % str(data[1])
                else:
                    # data[2] == derived results from point_processor
                    #   (0: parameters, 1: variables)
                    if var_count:
                        print "# Variables:"
                        print ", ".join(
                            "%s = %r" % (k, v)
                            for (k, v) in zip(var_names, data[1])
                        )

                    if data[2]:
                        print "# Data:"
                        print_table(
                            list(pad_list(data_names, "n/a", len(data[2]))),
                            data[2]
                        )

                    try:
                        output_columns = post_processor(
                            pars=data[0], vars=data[1], data=data[2], file=[]
                        )
                        print "# Output columns:"
                        print_table(
                            map(str, range(1, 1 + len(output_columns))),
                            output_columns
                        )
                    except Exception as e:
                        # in PARSPACE_MODE_FILE, this will surely raise an
                        #   error if file.NAME is used in out_columns
                        #   there is, however, no obvious way to get those
                        #   values
                        print "Error in 'out_columns':", e

                    if likelihood_code is not None or \
                       projection_count or symmetry_trafos:
                        # likelihood, projections or symmetries are extras
                        print "# Further information:"

                    if likelihood_code is not None:
                        try:
                            L = likelihood(data)
                            if min_likelihood is None:
                                print "likelihood =", L
                            else:
                                print "likelihood =", L,
                                # use "not a >= b" instead of "a < b" to
                                #   catch "nan"s and replicate the explorer
                                #   behavior
                                print "<" if not L >= min_likelihood else ">",
                                print min_likelihood

                            data[2].append(L)
                        except Exception as e:
                            print "Error in 'likelihood':", e
                            data[2].append(float("nan"))

                    def projection_eval(s):
                        """Evaluate projection formula."""
                        return formula_eval(
                            s,
                            pars=data[0],
                            vars=data[1],
                            data=data[2],
                            likelihood=likelihood(data)
                        )
                    for i, proj in enumerate(projections):
                        try:
                            (
                                x_formula, y_formula,
                                z_formula, filter_formula
                            ) = proj

                            x = projection_eval(x_formula)
                            y = projection_eval(y_formula)
                            z = projection_eval(z_formula)
                            if filter_formula:
                                if projection_eval(filter_formula):
                                    f_status = "(passed filter)"
                                else:
                                    f_status = "(failed filter)"
                            else:
                                f_status = "(no filter)"

                            print "projection_%i =" % i, (x, y, z), f_status
                        except Exception as e:
                            print "Error in projection %i: %s" % (i, str(e))

                    if symmetry_trafos:
                        print "# Transformed points:"
                        for i, sym in enumerate(symmetry_trafos):
                            try:
                                transformed_point = apply_symmetry(
                                    sym, point, par_names
                                )
                                print "%i:" % i, ", ".join(
                                    "%s = %r" % (k, v)
                                    for (k, v) in zip(
                                        par_names,
                                        transformed_point
                                    )
                                )
                            except FormulaError as e:
                                # only print the first part of the
                                #   FormulaError, since the second only
                                #   specifies the point (which we already
                                #   printed above)
                                print "%i:" % i, e.args[0]

                    # write calculated data into ".testdata" file for
                    #   possible later use
                    with DelayedKeyboardInterrupt(), \
                            open("%s.testdata" % base_name, "a") as f:
                        f.write("\t".join(map(out_repr, list_sum(data))))
                        f.write("\n")

                print
                print test_record_separator()
            except FormulaError as e:
                print e
            except EOFError as e:
                # see enclosing try/except block
                print
                exit_program()
            except Exception as e:
                # in case of an unexpected exception, i.e. not caught in
                #   process_item or similar, print out extended traceback
                traceback.print_exc()

    except (EOFError, KeyboardInterrupt) as e:
        # test mode has barely any setup, so just exiting is fine
        # also print one last \n for cosmetic purposes
        print
        exit_program()

elif mode == MODE_SCAN:
    # -------------------- #
    # |  [G.2] Scan Mode | #
    # -------------------- #
    if parspace_mode == PARSPACE_MODE_GRID:
        for r in par_ranges:
            if not r.is_finite:
                exit_program("Error: parameter range not finite: " + str(r))

        full_count = prod(map(len, par_ranges))
        # itertools.product does not allocate full product as intermediate,
        #   so this is safe to generate
        all_items = itertools.product(*par_ranges)
        # columns_mapping is a wrapper to the input columns
        #   (for parspace mode "grid", it's trivial)
        columns_mapping = lambda x: x
    elif parspace_mode == PARSPACE_MODE_SCATTER:
        full_count = scatter_point_count
        if all(p.is_finite for p in par_ranges):
            all_count = prod(map(len, par_ranges))
            if full_count > all_count:
                warn(
                    "Warning: scatter point count is higher than number",
                    "of possible points"
                )
        all_items = scatter_iterator(full_count, par_ranges)
        # it's also trivial for scatter scans
        columns_mapping = lambda x: x
    elif parspace_mode == PARSPACE_MODE_FILE:
        columns = []
        for i in range(par_count):
            columns.append(par_definitions[i])

        full_count = linecount(parspace_files)
        all_items = parspace_file_iterator(parspace_files)
        # columns_mapping is a wrapper to the input columns
        columns_mapping = lambda x: get_columns(columns, x)
    else:
        exit_program("Internal error: unknown parameter space mode")

    # read in resume position and skip ahead that many items
    try:
        with open(base_name + ".resume", "r") as f:
            done_count = int(f.read())
    except Exception:
        done_count = 0
    consume(done_count, all_items)

    warn(
        "# Parameter space mode:",
        dict(
            zip(parspace_mode_names.values(), parspace_mode_names.keys())
        )[parspace_mode]
    )
    warn("# Resuming at", done_count, "of", full_count)

    wait_for_user("start the scan")

    # number of valid data points already calculated
    #   this is not the same as the number of ALL data points calculated!
    found_count = linecount(base_name + ".data")

    # cache for overall calculation speed
    rates = []
    try:
        with open(base_name + ".speed") as f:
            for line in f:
                rates.append(float(line))
    except Exception:
        pass

    # each batch of data points has length <unit_length>
    batch = take(unit_length, all_items)
    progress_milestones = []
    try:
        draw_progress(done_count, full_count)
        while done_count < full_count:
            with ConfirmedExitOnInterrupt():
                # on user interrupt, batch has either been saved completely
                #   or can be started over without problems
                starttime = time.time()

                # explicitly unravel iterator provided by calculate_items
                # keep track of initial data before columns_mapping to have
                #   file data at hand in items mode "file"
                # NOTE: itertools.product returns a tuple, but everything
                #   else is a list so we have to convert "item" to list
                # in grid mode item is just the tuple of parameter values, in
                #   files mode it is the tuple of input file column values;
                #   so in principle pars[i] == file[i] in grid mode
                points_with_data = list(
                    calculate_items([
                        (columns_mapping(item), list(item)) for item in batch
                    ])
                )

                # delay interrupts until after saving of state, so once we
                #   calculate something it is saved in full and no repeats or
                #   missing points happen
                with DelayedKeyboardInterrupt(
                        "# Saving data to disk, please wait."
                        ):
                    # also to counteract possible errors in the formula used
                    #   in post_processor, we first build all the rows and
                    #   only then save them to the files; this way there are
                    #   no half-done batches
                    data_rows = []
                    excluded_rows = []
                    for point_and_info in points_with_data:
                        (point, old_data) = point_and_info
                        # put valid data in ".data" and invalid in
                        #   ".excluded-data"
                        if not isinstance(point, Exception):
                            with open(base_name + ".data", "a") as f:
#                                data_rows.append( point[2] )

                                data_rows.append(post_processor(
                                    pars=point[0],
                                    vars=point[1],
                                    data=point[2],
                                    file=old_data
                                ))

#                                data_rows.append(post_processor(
#                                    pars=point[0],
#                                    vars=point[1],
#                                    data=point[2],
#                                    file=old_data
#                                ))
                            found_count += 1
                        # Warning/WARNING/warning
                        # - Commented out not to have huge excluded data files
#                       else:
#                           with open(base_name + ".excluded-data", "a") as f:
#                               excluded_rows.append(
#                                   [ERROR_MARKER] +
#                                   map(out_repr, point.args[0]) +
#                                   [point.args[1]]
#                               )

                    with open(base_name + ".data", "a") as f:
                        for row in data_rows:
                            f.write("\t".join(row))
                            f.write("\n")

#Added by Ciara
                    with open("top_dir_/save.data", "a") as f:
                        for row in data_rows:
                            f.write("\t".join(row))
                            f.write("\n")


                    # Warning/WARNING/warning
                    # - Commented out not to have huge excluded data files

#                   with open(base_name + ".excluded-data", "a") as f:
#                       for row in excluded_rows:
#                           f.write("\t".join(row))
#                           f.write("\n")

                    done_count += len(batch)
                    with open(base_name + ".resume", "w") as f:
                        f.write(repr(done_count))

                    batch = take(unit_length, all_items)

                # cache calculated speed for later restarts of program
                # speed = mean rate of points/second over full time frame of
                #   speed statistics
                endtime = time.time()
                rate = unit_length / (endtime - starttime)
                rates.append(rate)
                with open(base_name + ".speed", "a") as f:
                    f.write("%f" % rate)
                    f.write("\n")

                rate = mean(rates)
                rateerror = std(rates)
                minrate = max(rate - rateerror, min(rates))

                # milestone == progress in %, rounded down to nearest 10%
                milestone = 10 * ((done_count * 10) // full_count)
                if milestone not in progress_milestones:
                    log_msg("# %i%% done" % milestone)
                    progress_milestones.append(milestone)

                draw_progress(
                    done_count, full_count,
                    "%d done, %d valid\nspeed: %f +- %f / s, ETA: %s to %s" % (
                        done_count,
                        found_count,
                        rate, rateerror,
                        progress_forecast(done_count, full_count, rate),
                        progress_forecast(done_count, full_count, minrate)
                    )
                )
    except Exception as e:
        # Errors during the scan abort it immediately to minimize damage
        # Aborts due to errors during the saving of errors can lead to
        #   inconsistency between resume value and number of points in the
        #   data files. However, in the case of errors in user supplied
        #   formulas, all data should considered tainted anyway and be redone
        #   (except in the case where file access rights changed during the
        #   loop and effected the error state, but that is not our fault)
        if cli_arguments.debug and not isinstance(e, FormulaError):
            traceback.print_exc()
        exit_program("Error during scan: " + str(e))

elif mode == MODE_MCMC:
    # -------------------- #
    # |  [G.3] MCMC Mode | #
    # -------------------- #

    # NOTE: MCMC does not support workers, but one can start MCMC on other
    #   computers independently without resorting to manager/worker
    #   architecture anyway

    chain_count = processing_concurrency()
    points_with_data = [None for i in range(chain_count)]

    lengths = [0 for i in range(chain_count)]
    iterations = [0 for i in range(chain_count)]
    warn("# Finding starting points")

    # find suitable starting points for all chains
    for i in range(chain_count):
        length = 0
        iteration = 0
        # first try to get some from previous calculations
        # (note that finding an invalid(or likelihood=0) point in the valid
        #   data should never happen and is an error)
        try:#Change by Ciara
            with open(base_name + ".chain.%i" % i) as f and open("top_dir_/save.chain.%i" i):
                last_line = ""
                for line in f:
                    length += 1
                    iteration += int(line.split("\t")[-1])
                    # last column is always stay count
                    last_line = line

                if not length:
                    continue

                full_data = map(float, last_line.split("\t"))
#Added by Ciara 
                print("full_data", full_data)
                if full_data[0] == ERROR_MARKER:
                    raise ValueError("invalid points in chain file")
                # inf == infinity will cause all columns after the first two
                #   groups to be put into group 3
                # thus we have to first manually drop the last two columns
                #   (likelihood and staycount)
                point_with_data = split_cols(
                    full_data[0:-2],
                    [par_count, var_count, float('inf')]
                )
#Added by Ciara
                print("point_with_data ", point_with_data)
                lengths[i] = length
                iterations[i] = iteration

                L = likelihood(point_with_data)
                if L == 0:
                    # either the likelihood function changed or something
                    #   went wrong
                    # either way we cannot continue with this chain
                    raise ValueError("chain %i has point with likelihood 0")

                points_with_data[i] = point_with_data
        except IOError as e:
            # file cannot be opened
            pass
        except ValueError as e:
            # can e.g. catch errors in int() of iteration reading, or float
            #   in full_data, or simply a chain with points with L == 0
            if cli_arguments.debug:
                traceback.print_exc()
            exit_program("Error: " + str(e))
        except KeyboardInterrupt:
            exit_program()

    # if some chains are missing starting points, find some
    missing = it_len(x for x in points_with_data if x is None)
    try:
        errors = 0
        while missing:
            # take some random candidates...
            candidates = [
                make_random_point(par_ranges) for i in range(missing)
            ]
            # ... and pass them to the process pool
            for point in calculate_items(candidates):
                if not isinstance(point, Exception):
                    L = likelihood(point)
                    if L > 0:
                        # we have a good point
                        #   -> put it into the first empty slot
                        for i in range(len(points_with_data)):
                            if points_with_data[i] is not None:
                                continue
                            points_with_data[i] = point
                            missing -= 1
                            break
                    elif L < 0:
                        raise ValueError(
                            "likelihood is negative for point: " +
                            repr(point[0])
                        )
                else:
                    errors += 1

                draw_progress(
                    chain_count - missing,
                    chain_count,
                    "Encountered errors: %d" % errors
                )
    except KeyboardInterrupt:
        # if the user decides to interrupt here, just exit
        print
        exit_program()
    except (FormulaError, ValueError) as e:
        if cli_arguments.debug and not isinstance(e, FormulaError):
            traceback.print_exc()
        exit_program("Error: " + str(e))

    warn("# -> Done")
    warn("# Resuming from", lengths, "to", unit_length)

    wait_for_user("start")
    chain_seeds = [random.random() for i in range(chain_count)]
#Ciara changed here:
    #with TemporaryDirectory() as chain_status_dir, ConfirmedExitOnInterrupt():
    with top_dir_/chain_status_dir as chain_status_dir, ConfirmedExitOnInterrupt():
        # MCMC mode does not use process_batch but rather its pool directly
        #   calling map_async instead of calculate_items or process_batch or
        #   process_item
        # (MCMC does not make use of workers)
        # ParameterRange is picklable, so no problem here
        result = process_batch.pool.map_async(
            make_chain,
            [
                (
                    os.path.abspath(base_name + ".chain.%i" % chain_id),
                    os.path.abspath(base_name + ".rejected.%i" % chain_id),
                    starting_point,
                    lengths[chain_id],
                    unit_length,
                    par_ranges,
                    iterations[chain_id],
                    likelihood_code,
                    heating_function_code,
                    os.path.join(
                        chain_status_dir, "chain.%i.status" % chain_id
                    ),
                    cli_arguments.debug,
                    chain_seeds[chain_id]
                )
                for
                (chain_id, starting_point) in enumerate(points_with_data)
            ]
        )

        deltat = 2.0  # sleeping time between progress update
        progress0 = sum(lengths)
        last_time = time.time()
        full_length = unit_length * chain_count
        progress_milestones = []
        draw_progress(progress0, full_length)

        rates = []
        try:
            with open(base_name + ".speed") as f:
                for line in f:
                    rates.append(float(line))
        except Exception:
            pass

        last_loop = False
        while True:
            time.sleep(deltat)
            lengths = [0] * chain_count
            iterations = [0] * chain_count
            warnings = []
            for i in range(chain_count):
                try:
                    chain_status_file = os.path.join(
                        chain_status_dir, "chain.%i.status" % i
                    )
                    with open(chain_status_file) as status_file:
                        length, iteration, stuck_reasons = pickle.load(
                            status_file
                        )
                        lengths[i] = length
                        iterations[i] = iteration
                    if sum(stuck_reasons.values()) > 100:
                        warnings.append(
                            ("chain %i has been stuck for %i iterations "
                             "due to %s") % (
                                i,
                                sum(stuck_reasons.values()),
                                ", ".join(
                                    "%s (%i)" % (r, i)
                                    for r, i in stuck_reasons.iteritems()
                                )
                            )
                        )
                except Exception as e:
                    pass

            progress1 = max(progress0, sum(lengths))

            # do not divide by deltat here, since the user could have kept us
            #   in the abort prompt for much longer
            rate = (progress1 - progress0) / (time.time() - last_time)
            progress0 = progress1
            last_time = time.time()
            rates.append(rate)
            with open(base_name + ".speed", "a") as f:
                f.write("%f" % rate)
                f.write("\n")

            rate = mean(rates)
            rateerror = std(rates)
            minrate = max(rate - rateerror, min(rates))

            # milestone == progress in %, rounded down to nearest 10%
            milestone = 10 * ((progress1 * 10) // full_length)
            if milestone not in progress_milestones:
                log_msg("# %i%% done" % milestone)
                progress_milestones.append(milestone)

            draw_progress(
                progress1, full_length,
                "chain lengths: %s\nspeed: %f +- %f / s, ETA: %s to %s%s" % (
                    ", ".join(map(str, lengths)),
                    rate, rateerror,
                    progress_forecast(progress1, full_length, rate),
                    progress_forecast(progress1, full_length, minrate),
                    "\n%s" % "\n".join(warnings) if warnings else ""
                )
            )

            if last_loop:
                break
            if result.ready():
                last_loop = True

    try:
        result.get()
    except Exception as e:
        if cli_arguments.debug and not isinstance(e, FormulaError):
            traceback.print_exc()
        exit_program("Error: " + str(e))

elif mode == MODE_OPTIMIZE:
    # ------------------------ #
    # |  [G.4] Optimize Mode | #
    # ------------------------ #
    range_is_iterable = [False] * par_count
    for i in range(par_count):
        r = par_ranges[i]
        range_is_iterable[i] = r.is_finite

    cache = dict()
    enable_cache = all(range_is_iterable[i] for i in range(par_count))

    population = set()
    pre_population = set()

    # ".population" contains likelihood as first column
    #   rest of the columns are the coordinates
    try:
        if os.path.exists(base_name + ".population"):
            with open(base_name + ".population") as f:
                for line in f:
                    cols = map(float, line.split("\t"))
                    x = list(cols[1:])

                    for i, r in enumerate(par_ranges):
                        if r.is_finite:
                            x[i] = r.truncate(x[i])
                    x = tuple(x)

                    if len(x) != par_count:
                        raise ValueError(
                            "malformed line in population file: " + line
                        )
                    L = cols[0]
                    population.add((x, L))
                    if enable_cache:
                        cache[x] = L
    except (IOError, ValueError) as e:
        warn("Discarding population file due to error: ", e)
        population = set()
        cache = dict()

    if population:
        warn("# Loaded population from file:", len(population), "points")

    # fill cache and pre_population with ".data" and ".excluded-data",
    #   but only if we need it
    if enable_cache or len(population) < unit_length:
        try:
            with ConfirmedExitOnInterrupt():
                for suffix in [".data", ".excluded-data"]:
                    if os.path.isfile(base_name + suffix) and \
                       os.path.getsize(base_name + suffix):
                        warn(
                            "# Filling cache and population candidates from",
                            "previous data:",
                            base_name + suffix
                        )
                        file_size = os.path.getsize(base_name + suffix)
                        read_size = 0
                        with open(base_name + suffix) as f:
                            for line in f:
                                read_size += len(line)
                                cols = map(maybefloat, line.split("\t"))
                                if cols[0] == ERROR_MARKER:
                                    # excluded lines begin with ERROR_MARKER
                                    #   and are only pars + error after that
                                    L = float("nan")
                                    x = tuple(cols[1:par_count + 1])
                                else:
                                    # inf == infinity will cause all columns
                                    #   after the first two groups to be put
                                    #   into group 3
                                    data = split_cols(
                                        cols,
                                        [par_count, var_count, float('inf')]
                                    )
                                    L = likelihood(data)
                                    x = cols[:par_count]
                                    for i, r in enumerate(par_ranges):
                                        if r.is_finite:
                                            x[i] = r.truncate(x[i])

                                    x = tuple(x)

                                if enable_cache:
                                    cache[x] = L
                                if not math.isnan(L):
                                    pre_population.add((x, L))

                                draw_progress(
                                    read_size, file_size,
                                    ("population candidates: %i; "
                                     "cached values: %i") % (
                                        len(pre_population),
                                        len(cache)
                                    )
                                )
        except Exception as e:
            if cli_arguments.debug and not isinstance(e, FormulaError):
                traceback.print_exc()
            exit_program("Error during loading of available data: " + str(e))

    try:
        if pre_population and len(population) < unit_length:
            answer = raw_input(
                "# You have already calculated points, but have no current "
                "full population - should we use those points to fill up to "
                "full? [Y/n]"
            )
            print "# Adding most optimal point found so far to population"
            population.add(max(pre_population, key=lambda x_and_L: x_and_L[1]))
            if answer.lower() != "n":
                if len(pre_population) >= unit_length - len(population):
                    answer = raw_input(
                        "# Should the remaining points be chosen at random "
                        "(Y) or from the most fit (n)? [Y/n]"
                    )
                    if answer.lower() != "n":
                        # choose randomly
                        population.update(
                            random.sample(
                                pre_population,
                                unit_length - len(population)
                            )
                        )
                    else:
                        # choose from the best
                        population.update(
                            sorted(
                                pre_population,
                                key=lambda x: -x[1]
                            )[:unit_length - len(population)]
                        )
                else:
                    # pre_population is so small that random and best choice
                    #   is equivalent
                    population.update(pre_population)

    except (KeyboardInterrupt, EOFError) as e:
        print
        exit_program()

    if not population:
        warn("# No previous population found")
    wait_for_user("start")

    # make up a population of unit_length points
    if len(population) < unit_length:
        warn("# Filling population to", unit_length, "points")

    try:
        with ConfirmedExitOnInterrupt():
            error_count = 0
            while len(population) < unit_length:
                # if we are lacking some, try new random candidates but only
                #   as many as necessary
                candidates = set()
                while len(candidates) < unit_length - len(population):
                    candidates.add(tuple(make_random_point(par_ranges)))
                candidates = sorted(candidates)

                # only save progress info over the last 1000 point bunches
                times = deque(
                    [time.time()],
                    maxlen=1000 * concurrent_processors
                )
                for data in calculate_items(candidates):
                    # make the handling of each point atomic w.r.t. interrupts
                    with DelayedKeyboardInterrupt():
                        times.append(time.time())
                        if not isinstance(data, Exception):
                            L = likelihood(data)
                            dressed_data = (tuple(data[0]), L)
                            # also append L to result columns
                            data[-1].append(L)

                            population.add(dressed_data)
                            if enable_cache:
                                cache[dressed_data[0]] = dressed_data[1]
                            suffix = ".data"
                        else:
                            data = [
                                [ERROR_MARKER],
                                list(data.args[0]),
                                [data.args[1]]
                            ]
                            suffix = ".excluded-data"
                            error_count += 1

                        with open(base_name + suffix, "a") as f:
                            f.write("\t".join(map(out_repr, list_sum(data))))
                            f.write("\n")

                        rate, rateerror, minrate = times_to_rate_info(
                            times,
                            concurrent_processors
                        )

                        draw_progress(
                            len(population), unit_length,
                            (
                                "speed: %f +- %f / s, ETA: %s to %s\n"
                                "errors: %i"
                            ) % (
                                rate, rateerror,
                                progress_forecast(
                                    len(population),
                                    unit_length,
                                    rate
                                ),
                                progress_forecast(
                                    len(population),
                                    unit_length,
                                    minrate
                                ),
                                error_count
                            )
                        )

    except Exception as e:
        if cli_arguments.debug and not isinstance(e, FormulaError):
            traceback.print_exc()
        exit_program("Error during loading: " + str(e))

    population = sorted(population)
    warn(
        "# Starting with population with fitness:",
        repr(max(x[1] for x in population))
    )

    waiting = 0
    # wait at least a number of iterations where the difference stays smaller
    #   than epsilon and only then consider the process to have converged
    try:
        with ConfirmedExitOnInterrupt():
            while waiting <= max_waiting_time:
                # ".population" contains likelihood as first column, rest of
                #   the columns are coordinates
                with DelayedKeyboardInterrupt(
                        "# Saving data to disk, please wait."
                        ):
                    with open(base_name + ".population", "w") as f:
                        for x in population:
                            f.write("\t".join(map(out_repr, (x[1], ) + x[0])))
                            f.write("\n")

                done = False
                while not done:
                    new_population = list()
                    for x, L in population:
                        triple = []
                        # choose random distinct a, b, c from current
                        #   population
                        while len(triple) < 3:
                            a = random.choice(population)[0]
                            if a != x and a not in triple:
                                triple.append(a)

                        y = []
                        for i, r in enumerate(par_ranges):
                            # three parent mating formula: xs = a + S (b - c)
                            # S being differential_weight
                            if range_is_iterable[i]:
                                # if the range is iterable, understand
                                #   formula as applying to indices into list
                                #   of values
                                xs_i = int(round(
                                    r.index(triple[0][i]) +
                                    differential_weight * (
                                        r.index(triple[1][i]) -
                                        r.index(triple[2][i])
                                    )
                                ))
                                xs_i = r.values[min(len(r) - 1, max(xs_i, 0))]
                            else:
                                # otherwise use values as is
                                xs_i = triple[0][i] + differential_weight * (
                                    triple[1][i] - triple[2][i]
                                )
                                # but make sure it does not leave the
                                #   parameter range
                                xs_i = r.pull_inside(xs_i)

                            # candidate: randomly choose coordinate value
                            #   from x or xs
                            y.append(
                                x[i]
                                if random.random() < crossover_probability
                                else xs_i
                            )
                        y = tuple(y)
                        new_population.append(y)

                    done = True
                # in the event of an interrupt, continue by forfeiting
                #   new_population and redoing the evolution step loop

                # population is form (x, f(x)), new_population still only y
                pairs = zip(new_population, population)
                work_items = []
                processed_items = []
                # if we have cached values, use those for the new candidates
                for pair in pairs:
                    if enable_cache and pair[0] in cache:
                        cached_value = cache[pair[0]]
                        if math.isnan(cached_value):
                            # if cached likelihood is "nan", i.e. calculation
                            #   yielded an error, we have no good candidates
                            #   -> do no replacement
                            processed_items.append((pair[1], pair[1]))
                        else:
                            processed_items.append((
                                (pair[0], cached_value),
                                pair[1]
                            ))
                    else:
                        work_items.append(pair)

                count = 0
                # only save progress info over the last 1000 point bunches
                times = deque(
                    [time.time()],
                    maxlen=1000 * concurrent_processors
                )
                # process those points that have no cached values
                for new_data, old_point in calculate_items(work_items):
                    # delay interrupts, such that each point is handled
                    #   atomically
                    with DelayedKeyboardInterrupt():
                        times.append(time.time())
                        count += 1
                        if not isinstance(new_data, Exception):
                            L = likelihood(new_data)
                            if enable_cache:
                                cache[tuple(new_data[0])] = L

                            new_point = (tuple(new_data[0]), L)
                            processed_items.append((new_point, old_point))

                            # also append L to result columns
                            new_data[-1].append(L)
                            suffix = ".data"
                        else:
                            if enable_cache:
                                cache[tuple(new_data[0])] = float("nan")
                            processed_items.append((old_point, old_point))

                            new_data = [
                                [ERROR_MARKER],
                                list(new_data.args[0]),
                                [new_data.args[1]]
                            ]
                            suffix = ".excluded-data"

                        with open(base_name + suffix, "a") as f:
                            f.write("\t".join(
                                map(out_repr, list_sum(new_data))
                            ))
                            f.write("\n")

                        rate, rateerror, minrate = times_to_rate_info(
                            times,
                            concurrent_processors
                        )

                        draw_progress(
                            count, len(work_items),
                            "speed: %f +- %f / s, ETA: %s to %s" % (
                                rate, rateerror,
                                progress_forecast(
                                    count,
                                    len(work_items),
                                    rate
                                ),
                                progress_forecast(
                                    count,
                                    len(work_items),
                                    minrate
                                )
                            )
                        )

                with DelayedKeyboardInterrupt(
                        "# Making new population, please wait."
                        ):
                    merged_population = []
                    for y, x in processed_items:
                        # only take y (new) instead of x if it's better
                        # NOTE: this also falls back to x if y[1] is "nan",
                        #   because a > nan is always False for all a
                        merged_population.append(y if y[1] > x[1] else x)

                    max_old = max(x[-1] for x in population)
                    max_new = max(x[-1] for x in merged_population)
                    delta = abs(max_old - max_new)
                    delta_max = epsilon + epsilon_relative * abs(max_old)
                    if delta <= delta_max:
                        waiting += 1
                    else:
                        waiting = 0

                    # if |merged_population| < |population|, fill up with the
                    #   best points from both new and old
                    #   i.e. with the data contained in processed_items
                    merged_population = set(merged_population)
                    missing = len(population) - len(merged_population)
                    filled_up = 0
                    if filled_up < missing:
                        candidates = []
                        for y, x in processed_items:
                            candidates.append(x)
                            candidates.append(y)

                        for x in sorted(candidates, key=lambda x: x[1]):
                            if x not in merged_population:
                                merged_population.add(x)
                                filled_up += 1
                                if filled_up >= missing:
                                    break
                        else:
                            # for loop was not broken out of
                            #   => there are still missing items
                            raise RuntimeError(
                                "Internal error: "
                                "population cannot be replenished"
                            )

                    population = sorted(merged_population)
                    if waiting > 1:
                        stability_msg = "(stable for %i iterations)" % (
                            waiting - 1
                        )
                    else:
                        stability_msg = ""

                    warn(
                        "# New population with fitness:",
                        repr(max([x[1] for x in population])),
                        stability_msg
                    )
                    # when interrupt happens, the block will still finish and
                    #   everything will be fine

                with DelayedKeyboardInterrupt(), \
                        open(base_name + ".optimum", "w") as f:
                    optimum = max(population, key=lambda x_and_L: x_and_L[1])
                    f.write("\t".join(
                        map(out_repr, (optimum[1], ) + optimum[0])
                    ))
                    f.write("\n")

    except Exception as e:
        if cli_arguments.debug and not isinstance(e, FormulaError):
            traceback.print_exc()
        exit_program("Error during evolution process: " + str(e))

    exit_program(
        "# Maximal waiting time exceeded => algorithm converged",
        successful=True
    )

elif mode == MODE_EXPLORER:
    # ------------------------ #
    # |  [G.5] Explorer Mode | #
    # ------------------------ #

    for i, r in enumerate(par_ranges):
        if not r.is_finite:
            exit_program("Error: parameter range not finite: " + str(r))

    def truncate_pars(pars):
        """Truncate all parameter values to their respective ranges."""
        return tuple(r.truncate(pars[i]) for i, r in enumerate(par_ranges))

    full_length = prod(map(len, par_ranges))

    print "# In total %i possible points" % full_length

    # we save as few values as necessary in memory on each point
    # necessary fields: likelihood, maybe_boundary, projection coordinates,
    #   projection filters
    # (maybe_boundary is initially True and is set to False once all
    #   neighbors have been calculated)
    saved_fields = ["likelihood", "maybe_boundary"]
    # projections work by grouping points by two coordinates and then
    #   maximizing the third coordinate for every such x-y-equivalent group
    # z coordinate is either explicitly taken from config or likelihood per
    #   default
    # if the config option "projection_#" is given directly, only parameters
    #   and variables can be used in x, y coordinates
    projection_operators = []

    # decide when coordinate or filter result has to be cached during dressing
    #   this must never yield false negatives (uncached complicated formulas)
    depends_on_vars_or_data = lambda s: "data" in s or "vars" in s
    is_complicated = depends_on_vars_or_data
    for proj in projections:
        xyzF_formulas = list(proj)

        for i, f in enumerate(xyzF_formulas):
            if f and is_complicated(f):
                if f not in saved_fields:
                    saved_fields.append(f)
                xyzF_formulas[i] = "fields[%i]" % saved_fields.index(f)

        projection_operators.append(
            "(%s, %s, %s)%s" % (
                xyzF_formulas[0],
                xyzF_formulas[1],
                xyzF_formulas[2],
                (" if (%s) else None" % xyzF_formulas[3])
                if xyzF_formulas[3] else ""
            )
        )
    saved_fields = tuple(saved_fields)

    # dressed data = likelihood + maybe_boundary + any fields needed for
    #   projections

    # maybe_boundary is set to False, if all neighbors have been calculated
    # NOTE: out of the three sqlite, list, array for data storage, array has
    #   the best performance while also being the most lightweight in memory
    #   usage
    def dress_data(data, saved_fields=saved_fields, likelihood=likelihood):
        """Return derived point values necessary for exploration.

        Return array-type object that contains the (hopefully) bare minimum
            of data needed for the algorithm.
        """
        L = likelihood(data)
        return array.array(
            "f",
            [L, True] +
            [
                formula_eval(
                    field,
                    pars=data[0],
                    vars=data[1],
                    data=data[2],
                    likelihood=L
                ) for field in saved_fields[2::]
            ]
        )

    # nodes gives index -> dressed data relation
    # NOTE: this is the only info on all points that is kept in memory, full
    #   data is ONLY saved in ".data" file
    # Thus expanding the set of projections might force a regeneration of
    #   this data set
    nodes = dict()

    input_files = parspace_files + [
        base_name + ".data", base_name + ".excluded-data"
    ]
    if any(map(os.path.isfile, input_files)):
        warn("# Importing data")

    # interesting subset mechanics:
    # subnodes contains every point index that was in the last work load or
    #   to-be-explored boundary set, i.e. points that have been touched very
    #   recently
    # if that file exists and the user agrees, interestingsubset is filled
    #   with its contents and all neighbors of those points
    # following that, from the whole set of calculated points only those are
    #   loaded into memory that are in that interesting subset or are erroneous
    # those selected points (if selection happens) are then written to the
    #   subdata file during loading to facilitate faster loading next time
    interestingsubset = set()

    if cli_arguments.experimental and sys.stdin.isatty() and \
       os.path.isfile(base_name + ".subnodes"):
        try:
            answer = raw_input(
                "# Do you want to skip points that are not "
                " connected to the last work set? [N/y]"
            )
            if answer != "y":
                raise TypeError("User did not type yes")
            interestingsubset = set()
            with open(base_name + ".subnodes", "r") as f:
                for line in f:
                    point_id = int(line)
                    interestingsubset.add(point_id)
                    interestingsubset.update(neighbors(point_id, par_ranges))

        except Exception:
            interestingsubset = set()

    full_unlikely_count = 0
    for filename in input_files:
        try:
            if not os.path.isfile(filename):
                continue
            with open(filename) as f:
                errors = 0
                error_msgs = []
                count = 0
                valid = 0
                duplicates = 0
                filtered = 0
                unlikely = 0
                full_size = os.path.getsize(filename)
                read_size = 0
                for line in f:
                    count += 1
                    read_size += len(line)

                    try:
                        cols = map(maybefloat, line.strip().split("\t"))
                        if cols[0] == ERROR_MARKER:
                            # excluded lines start with a column containing
                            #   ERROR_MARKER and are pars + error after that
                            data = split_cols(cols[1:], [par_count, 1])
                            pars = truncate_pars(data[0])
                            data = Exception(
                                pars,
                                data[1] if data[1] else "Unknown error"
                            )
                        else:
                            # inf = infinity will cause all columns after the
                            #   first two groups to be put into group 3
                            data = split_cols(
                                cols,
                                [par_count, var_count, float('inf')]
                            )
                            pars = truncate_pars(data[0])
                            # re-check bounds just to be safe
                            try:
                                check_bounds(bounds, pars, data[1], data[2])
                            except Exception as e:
                                data = Exception(pars, str(e))

                        if not isinstance(data, Exception) and \
                           not loading_filter(data):
                            filtered += 1
                            continue

                        index = point2index(pars, par_ranges)

                        # as mentioned only keep track of interesting (if
                        #   applicable) and invalid points
                        if index not in nodes:
                            if not interestingsubset or \
                               index in interestingsubset or \
                               isinstance(data, Exception):
                                # add point to nodes table if there's no
                                #   interestingsubset or
                                #   if it's in it or
                                #   if it's an Exception anyway

                                # only dress_data if it's actually valid
                                if not isinstance(data, Exception):
                                    dressed_data = dress_data(data)
                                    if not dressed_data[0] >= min_likelihood:
                                        # use "not a >= b" instead of
                                        #   "a < b" to catch "nan"s
                                        unlikely += 1
                                        full_unlikely_count += 1
                                    nodes[index] = dressed_data
                                else:
                                    nodes[index] = None

                            # write the interesting subset to ".subdata"
                            if interestingsubset:
                                with DelayedKeyboardInterrupt(), \
                                        open(base_name + ".subdata", "a") as f:
                                    f.write("\t".join(map(out_repr, cols)))
                                    f.write("\n")

                        else:
                            duplicates += 1

                        if not isinstance(data, Exception):
                            valid += 1
                    except ValueError as e:
                        # this usually catches errors from list.index in
                        #   point2index, i.e. points not lying on grid
                        if cli_arguments.debug:
                            error_msgs.append(e)
                        errors += 1

                    draw_progress(
                        read_size, full_size,
                        ("%i errors, %i imported, %i import valid, "
                         "%i duplicates, %i filtered, %i unlikely\n%s%s") % (
                            errors,
                            len(nodes),
                            valid,
                            duplicates,
                            filtered,
                            unlikely,
                            "file: %s" % abbreviate_path(
                                filename,
                                width=terminal_size()[0] - 6
                            ),
                            ("\n" + "\n".join(map(str, error_msgs[-10:])))
                            if cli_arguments.debug and error_msgs else ""
                        )
                    )
        except Exception as e:
            # if we find an exception here, it must have come from errors in
            #   the files or errors in the formulas
            # either way the user must first clear this up before we can
            #   continue
            if cli_arguments.debug and not isinstance(e, FormulaError):
                traceback.print_exc()
            exit_program("Error during loading: " + str(e))
        except KeyboardInterrupt:
            # if user interrupts, then just exit
            print
            exit_program()

    del interestingsubset

    # see if there is previously begun work
    #   (in the form of a file with point indices as lines)
    try:
        work = set()
        if os.path.isfile(base_name + ".work"):
            with open(base_name + ".work") as f:
                for line in f:
                    work.add(int(line))

        if work:
            warn(
                "# Loaded previously determined work load:",
                len(work),
                "points"
            )
    except Exception as e:
        warn("Error loading work load: ", e)
        work = set()
    except BaseException:
        # left-over exception types are from sys.exit and ^C
        exit_program()

    wait_for_user("start")

    # if there are no nodes and there is no work, try to find some randomly
    if not work and len(nodes) - full_unlikely_count < unit_length:
        warn(
            "# Filling data set to", unit_length,
            "points (with above minimal likelihood)"
        )
    try:
        with ConfirmedExitOnInterrupt():
            error_count = 0
            while not work and len(nodes) - full_unlikely_count < unit_length:
                # try this until we have enough, because, even though we always
                #   calculate enough points to fill up, some might turn out to
                #   be invalid
                candidates = set()
                while len(candidates) < \
                        unit_length - len(nodes) + full_unlikely_count:
                    candidates.add(tuple(make_random_point(par_ranges)))
                candidates = sorted(candidates)
                candidates = [
                    (p, point2index(p, par_ranges)) for p in candidates
                ]
                candidates = [
                    (p, idx) for (p, idx) in candidates
                    if idx not in nodes
                ]

                # only save progress info over the last 1000 point bunches
                times = deque(
                    [time.time()],
                    maxlen=1000 * concurrent_processors
                )
                for data, index in calculate_items(candidates):
                    # make the handling of each point atomic w.r.t. interrupts
                    with DelayedKeyboardInterrupt():
                        times.append(time.time())
                        if not isinstance(data, Exception):
                            # point -> tuple of parameters
                            point = data[0]

                            L = likelihood(data)
                            if not L >= min_likelihood:
                                # use "not a >= b" instead of "a < b" to
                                #   catch "nan"s
                                continue
                                # full_unlikely_count += 1
                            # also append likelihood to data values
                            data[2].append(L)

                            index = point2index(point, par_ranges)
                            nodes[index] = dress_data(data)
                            del index
                            suffix = ".data"
                        else:
                            data = [
                                [ERROR_MARKER],
                                list(data.args[0]),
                                [data.args[1]]
                            ]
                            suffix = ".excluded-data"
                            error_count += 1

                        with open(base_name + suffix, "a") as f:
                            f.write("\t".join(map(out_repr, list_sum(data))))
                            f.write("\n")

                        rate, rateerror, minrate = times_to_rate_info(
                            times,
                            concurrent_processors
                        )

                        draw_progress(
                            len(nodes) - full_unlikely_count, unit_length,
                            (
                                "speed: %f +- %f / s, ETA: %s to %s\n"
                                "errors: %i"
                            ) % (
                                rate, rateerror,
                                progress_forecast(
                                    len(nodes) - full_unlikely_count,
                                    unit_length,
                                    rate
                                ),
                                progress_forecast(
                                    len(nodes) - full_unlikely_count,
                                    unit_length,
                                    minrate
                                ),
                                error_count
                            )
                        )
    except Exception as e:
        if cli_arguments.debug and not isinstance(e, FormulaError):
            traceback.print_exc()
        exit_program("Error: " + str(e))

    warn(
        "# Exploring data space starting from in-total",
        "%i points, %i being invalid" % (
            len(nodes),
            len(filter(lambda v: v is None, nodes.values()))
        )
    )

    # algorithm state:
    # 1: do first line of projections  (explore new points emanating from
    #   best point for each x, y pair)
    # 2: do second line of projections (explore new points emanating from
    #   best boundary point for each x, y pair)
    # 3: explore new points emanating from boundary points, i.e. valid points
    #   with missing neighbors
    # we start in the lowest possible state
    algorithm_state = min(available_states)
    algorithm_state_names = [
        "one (basic)", "two (extended)", "three (unprojected)"
    ]
    # after each run where no points are found, it resets back to the minimal
    #   allowed state

    boundary = set()
    tobeexplored_boundary = set()
    work_items = set()
    num_explored = 0
    first_iteration = True
    prev_algorithm_state = None
    # used to track how often in a row state one did not yield new points
    state_one_completeness = 0

    # the rest of the algorithm assumes that par_ranges are lists
    par_ranges = map(list, par_ranges)

    # this loop only exits if there is nothing more to calculate and the
    #   algorithm is in the maximal allowed state
    while True:
        try:
            with ConfirmedExitOnInterrupt():
                # on user interrupt:
                #   restart the whole loop again
                #   this is not overly problematic, since the previous work
                #   load logic will skip most of the work

                # expect KeyboardInterrupts
                #   interrupts during discovery phase will restart it
                #   work load determination is atomic, so once discovery is
                #   done the work load is available until
                #       been it's done
                # each point is handled atomically, so there are no
                #   half-done points
                # all other errors are considered fatal, since they most
                #   likely come from user supplied formulas
                #   or errors in file input or output
                if not projection_operators:
                    algorithm_state = 3

                starttime = time.time()
                extrapolated_items = set()
                projection_differences = dict()
                print "%s\n%s" % (horizontal_line(), horizontal_line())
                warn(
                    "# Starting iteration at",
                    time.strftime("%d.%m.%Y %H:%M"),
                    "in state",
                    algorithm_state_names[algorithm_state - 1]
                )
                print "# Current scan:", os.path.abspath(chosen_scan_setup)

                if not work and first_iteration:
                    # only show if no unfinished work load, but don't mess up
                    #   the code below
                    warn(
                        "# Projection algorithm state:",
                        algorithm_state_names[algorithm_state - 1]
                    )

                num_nodes = []
                # find points whose neighbors should be calculated
                if work:
                    # previous work load not finished
                    #   -> no need to find anything
                    pass
                elif algorithm_state == 3:  # and not work
                    # START OF STATE 3
                    # 2 or higher -> calculate boundary and use the
                    #   appropriate subset
                    # only load/calculate full boundary if we know it is
                    #   needed for determining new work items

                    # if there is an already calculated boundary set for
                    #   exactly our situation, load it!
                    # sha512 should hopefully be enough to minimize hash
                    #   collision problems
                    m = hashlib.sha512()
                    m.update(repr(min_likelihood))
                    for k in sorted(nodes.keys()):
                        m.update(repr(k))
                    boundary_cache_checksum = m.hexdigest()

                    boundary_from_cache = False
                    try:
                        with open(base_name + ".boundary", "r") as f:
                            checksum = f.readline().strip()
                            if checksum != boundary_cache_checksum:
                                raise ValueError("Wrong boundary cache")
                            warn("# Loading boundary from cache")
                            boundary = []
                            for line in f:
                                boundary.append(int(line))
                            boundary_from_cache = True
                    except Exception as e:
                        boundary = []

                    # only consider the following points for further
                    #   exploration:
                    # * ones that were previously on the boundary, but have
                    #     not been explored yet
                    # * ones that are new additions and are valid data
                    #     points
                    if not boundary_from_cache:
                        warn("# Calculating boundary (this can take a while)")
                        if not boundary:
                            possible_boundary = set(nodes.keys())
                        else:
                            possible_boundary = set(boundary)

                        # tobeexplored_boundary and work_items are not needed
                        #   anymore after this
                        # -> get rid of them immediately (before calculating
                        #   the boundary) and free memory
                        if tobeexplored_boundary:
                            possible_boundary = set(possible_boundary) - \
                                set(tobeexplored_boundary)
                            tobeexplored_boundary = set()
                        if work_items:
                            possible_boundary |= set(
                                index
                                for (UNUSED, index) in work_items
                                if nodes[index] is not None
                            )
                            work_items = set()
                        gc.collect()

                        boundary = [
                            index for index in possible_boundary
                            if is_boundary(index, nodes, par_ranges)
                        ]
                        boundary.sort()
                        # possible_boundary not needed anymore
                        possible_boundary = None
                        gc.collect()

                        # file will be cleaned up before the delay is dropped
                        with DelayedKeyboardInterrupt(), \
                                open(base_name + ".boundary", "w") as f:
                            f.write(boundary_cache_checksum)
                            f.write("\n")
                            for index in boundary:
                                f.write(repr(index))
                                f.write("\n")

                    # break up the boundary in steps given by
                    #   likelihood_steps (pre-sorted in descending order)
                    # and calculate <unit_length> many from the highest
                    #   likelihood bin starting with  the most likely ones
                    likelihood_bins = [[] for step in likelihood_steps]
                    for index in boundary:
                        L = nodes[index][0]
                        bin_no = first(
                            i for i in range(len(likelihood_steps))
                            if L >= likelihood_steps[i]
                        )
                        likelihood_bins[bin_no].append(index)

                    warn(
                        "# Size of boundary set:", len(boundary)
                    )
                    if len(likelihood_bins) > 1:
                        warn(
                            "# Likelihood bin occupancy:",
                            ", ".join(repr(len(b)) for b in likelihood_bins)
                        )
                    assert len(boundary) == sum(map(len, likelihood_bins))
                    if not boundary:
                        exit_program(
                            "# No boundary points -> nothing to do",
                            successful=True
                        )

                    (source_i, source_bin) = first(
                        x for x in enumerate(likelihood_bins) if x[1]
                    )
                    if len(likelihood_steps) > 1:
                        warn(
                            "# Now exploring boundary down to likelihood =",
                            likelihood_steps[source_i]
                        )
                    tobeexplored_boundary = set(take(
                        min(unit_length, len(source_bin)),
                        sorted(source_bin, key=lambda index: -nodes[index][0])
                    ))
                    num_nodes.append(len(tobeexplored_boundary))

                    # unload items from <nodes> if their likelihood falls
                    #   into an older (i.e. higher) bin than what is handled
                    #   at the moment
                    if source_i >= min_purge_age:
                        min_purge_likelihood = likelihood_steps[
                            source_i - min_purge_age
                        ]
                        warn(
                            "# Unloading all points with likelihood",
                            "higher than",
                            min_purge_likelihood,
                            "from memory"
                        )
                        for index in nodes.keys():
                            if nodes[index] is not None and \
                               nodes[index][0] >= min_purge_likelihood:
                                del nodes[index]
                        gc.collect()

                    # END OF COMPLETENESS STATE 2+
                else:  # algorithm_state in [1, 2] and not work
                    # START OF COMPLETENESS STATE 0/1
                    # extrapolated_items cleared beforehand

                    if prev_algorithm_state == algorithm_state \
                       and prev_algorithm_state == 1 and tobeexplored_boundary:
                        # if last state was one, then tobeexplored_boundary
                        #   is the set of projected points with z=z_max
                        #   thus if we are still doing a state one
                        #   iteration, the new z_max points are either one
                        #   of those or one the newly calculated ones
                        # NOTE: if tobeexplored_boundary is empty, the last
                        #   iteration was a resuming one instead
                        #   so we cannot use this (non-existent) information
                        new_or_best_items = list(itertools.chain(
                            tobeexplored_boundary,
                            (index for (point, index) in work_items)
                        ))
                        all_items = lambda items=new_or_best_items: (
                            (index, nodes[index]) for index in items
                        )
                        # this is already stored inside the all_items
                        #   definition, so we can delete it
                        del new_or_best_items
                        # tobeexplored_boundary and work_items are disjoint
                        #   as elements of work_items were not in nodes
                        #   beforehand and those on the boundary are
                        base_count = len(nodes) - len(work_items) - \
                            len(tobeexplored_boundary)
                    else:
                        # if any of the relevant iterations were state two,
                        #   either the tobeexplored_boundary was not the
                        #   z_max points or the boundary set might have
                        #   changed thereby excluding the former z_max
                        #   points from the loop
                        # either way we cannot take a shortcut
                        all_items = nodes.iteritems  # no () at the end!
                        base_count = 0

                    tobeexplored_boundary = set()
                    work_items = set()
                    gc.collect()

                    print horizontal_line()

                    for proj_id, projection in enumerate(projection_operators):
                        # saves association (x, y) -> (max_z, list of items)
                        plane = dict()

                        items = all_items()
                        count = base_count
                        items_count = len(nodes)

                        xyzF_formulas = list(projections[proj_id])
                        if xyzF_formulas[3]:
                            xyzF_formulas[3] = "if " + xyzF_formulas[3]
                        else:
                            xyzF_formulas[3] = ""
                        print "# Projection:", \
                            "{}, {}, {} {}".format(*xyzF_formulas)

                        for (index, fields) in items:
                            count += 1
                            draw_progress(
                                count,
                                items_count
                            )

                            if fields is None:
                                # invalid point -> ignored
                                continue

                            if not fields[0] >= min_likelihood:
                                # use "not a >= b" instead of "a < b" to
                                #   catch "nan"s
                                # unlikely point -> ignored
                                continue

                            if algorithm_state > 1 and \
                               not is_boundary(index, nodes, par_ranges):
                                continue

                            xyz = formula_eval(
                                projection,
                                pars=index2point(index, par_ranges),
                                fields=fields,
                                likelihood=fields[0]
                            )
                            if xyz is None:  # did not pass projection filter
                                continue

                            # convert xy to str so reading in old projection
                            #   file find the same numbers
                            xy, z = (repr(xyz[0]), repr(xyz[1])), xyz[2]

                            # projected x,y point not yet seen -> set it up
                            if xy not in plane:
                                plane[xy] = [z, [index]]

                            if plane[xy][0] < z:
                                plane[xy][0] = z
                                plane[xy][1] = [index]
                            elif plane[xy][0] == z:
                                plane[xy][1].append(index)

                        # calculate difference between old and new
                        #   projection plane and save the new plane to file
                        if algorithm_state == 1:
                            old_plane = dict()
                            projection_file = base_name + ".projection.%i" % \
                                proj_id
                            if os.path.isfile(projection_file):
                                with open(projection_file) as f:
                                    for line in f:
                                        x, y, z = pad_list(
                                            line.split("\t"), "0", 3
                                        )
                                        z = float(z)  # leave x, y as str
                                        old_plane[(x, y)] = [z]

                            # entries of difference: count of new points,
                            #   changes in known points
                            projection_differences[proj_id] = diff = [0, []]
                            for (k, v) in plane.iteritems():
                                if k not in old_plane:
                                    diff[0] += 1
                                else:
                                    v_old = old_plane[k][0]
                                    # also mind possibly lossy str conversion
                                    #   of floats here
                                    v_new = float(repr(v[0]))
                                    diff[1].append(v_new - v_old)

                            # if the new plane is worse than the old one,
                            #   something must have gone wrong
                            # -> invalidate info on difference (new points
                            #   count might still be accurate though)
                            if any(dz < 0 for dz in diff[1]):
                                diff[1] = [float("nan")]

                            del diff  # free up name, does not delete object

                            with DelayedKeyboardInterrupt(), \
                                    open(projection_file, "w") as f:
                                # save projection to file in x-y-z format
                                for (k, v) in plane.iteritems():
                                    # v is [z, [...indices...]]
                                    f.write(
                                        "\t".join(
                                            map(out_repr, [k[0], k[1], v[0]])
                                        )
                                    )
                                    f.write("\n")

                        old_count = len(tobeexplored_boundary)
                        num_nodes.append(len(plane))
                        for (zmax, subspace) in plane.values():
                            tobeexplored_boundary.update(subspace)

                        print "# -> found %i points (%i not yet in list)" % (
                            len(plane),
                            len(tobeexplored_boundary) - old_count
                        ), "maximizing z in projected plane"

                        # interpolation and extrapolation only makes sense
                        #   in state one
                        # in state 2, for each x,y the lists of points
                        #   with these projection coordinates are not ordered
                        # thus it makes no sense to single out the first point
                        if algorithm_state == 1 and \
                           proj_id in extrapolated_projections:
                            # find projection grid values, i.e. discrete
                            #   coordinates
                            # Note: the keys are strings, so we have to convert
                            #   them first to make sensible comparison
                            xs = sorted(
                                set(x for (x, y) in plane.keys()),
                                key=float
                            )
                            ys = sorted(
                                set(y for (x, y) in plane.keys()),
                                key=float
                            )

                            # for coordinate interpretation see below
                            # directions (some may be covered multiple times
                            #   with different end points):
                            # \.|./
                            # .\|/.
                            # --O--
                            # ./|\.
                            # /.|.\
                            # NOTE longer inter- and extrapolation
                            #   directions disabled to proving unnecessary
                            #   in tests
                            interpolation_pairs = [
                                ((-1, 0), (1, 0)),
                                ((0, -1), (0, 1)),
                                ((-1, -1), (1, 1)),
                                ((-1, 1), (1, -1)),

                                # ((-2, 0), (2, 0)), ((0, -2), (0, 2)),
                                # ((-2, -2), (2, 2)), ((-2, 2), (2, -2)),
                                # ((-1, 2), (1, -2)), ((-1, -2), (1, 2)),
                                # ((-2, -1), (2, 1)), ((-2, 1), (2, -1))
                            ]
                            extrapolation_pairs = [
                                ((1, 0), (0, 0)),
                                ((0, 1), (0, 0)),
                                ((-1, 0), (0, 0)),
                                ((0, -1), (0, 0)),

                                # ((-2, 2), (-1, 1)), ((0, 2), (0, 1)),
                                # ((2, 2), (1, 1)), ((2, 0), (1, 0)),
                                # ((-2, -2), (-1, -1)), ((0, -2), (0, -1)),
                                # ((-2, -2), (-1, -1)), ((-2, 0), (-1, 0)),
                                # ((2, 0), (0, 0)), ((0, 2), (0, 0)),
                                # ((-2, 0), (0, 0)), ((0, -2), (0, 0))
                            ]

                            # find all used relative x and y coordinates
                            (smoothen_range_x, smoothen_range_y) = map(
                                lambda x: sorted(set(x)),
                                [
                                    list_sum(
                                        [
                                            [pair[0][i], pair[1][i]]
                                            for pair in interpolation_pairs
                                        ] + [
                                            [pair[0][i], pair[1][i]]
                                            for pair in extrapolation_pairs
                                        ]
                                    ) for i in range(2)
                                ]
                            )

                            count = 0
                            old_count = len(extrapolated_items)

                            print "# Inter- and extrapolation:"
                            for xy in plane:
                                i = xs.index(xy[0])
                                j = ys.index(xy[1])

                                # build array of neighboring indices for
                                #   the current point
                                # Example:
                                # ( 0, 0) current point
                                # (-2, 1) labels the point having an x-value
                                #   2 steps lower and an y-value 1 step higher
                                neighborhood = dict()
                                for a in smoothen_range_x:
                                    for b in smoothen_range_y:
                                        if 0 <= i + a < len(xs) and \
                                           0 <= j + b < len(ys):
                                            (xp, yp) = (xs[i + a], ys[j + b])
                                            if (xp, yp) in plane:
                                                neighborhood[(a, b)] = \
                                                    plane[(xp, yp)][1]
                                                # foo[1] is list of indices

                                # (direction pairs, function)
                                inter_and_extrapolation = (
                                    (
                                        interpolation_pairs,
                                        interpolate_indices
                                    ),
                                    (
                                        extrapolation_pairs,
                                        extrapolate_indices
                                    )
                                )
                                for pairs, func in inter_and_extrapolation:
                                    for pair in pairs:
                                        for (i1, i2) in itertools.product(
                                            neighborhood.get(pair[0], []),
                                            neighborhood.get(pair[1], [])
                                        ):
                                            extrapolated_items.update(
                                                index
                                                for index in
                                                func(i1, i2, par_ranges)
                                                if index not in nodes
                                            )
                                count += 1

                                draw_progress(
                                    count,
                                    len(plane)
                                )
                            del i, j, xs, ys

                            print ("# -> found %i additional points "
                                   "interpolated/extrapolated from "
                                   "projected points") % (
                                len(extrapolated_items) - old_count
                            )

                        if algorithm_state == 1 and symmetry_trafos:
                            print "# Symmetries:"
                            # makes only sense if each subspace is sorted
                            #   (only in state 1)
                            old_count = len(extrapolated_items)
                            count = 0
                            error_count = 0
                            for sym in symmetry_trafos:
                                for xy, stuff in plane.iteritems():
                                    # stuff == [z value, list of points]
                                    point = list(
                                        index2point(stuff[1][0], par_ranges)
                                    )
                                    try:
                                        point = apply_symmetry(
                                            sym, point, par_names
                                        )
                                        idx = point2index(point, par_ranges)
                                        if idx not in nodes:
                                            extrapolated_items.add(idx)

                                    except Exception as e:
                                        error_count += 1

                                    count += 1

                                    draw_progress(
                                        count,
                                        symmetry_count * len(plane)
                                    )
                            print "# -> found", \
                                len(extrapolated_items) - old_count, \
                                "additional points from supposed symmetries"

                        del plane
                        gc.collect()
                        print horizontal_line()
                    # END OF COMPLETENESS STATE 0/1
                    warn(
                        "# Sizes of projection planes:",
                        ", ".join(map(repr, num_nodes))
                    )

                    if algorithm_state == 1:
                        projection_differences = [
                            projection_differences[i]
                            for i in range(projection_count)
                        ]
                        warn("# Change in projections:")
                        max_ = lambda x: max(x) if x else float("nan")
                        points_added = [
                            diff[0] for diff in projection_differences
                        ]
                        norm1 = [
                            sum(abs(dz) for dz in diff[1])
                            for diff in projection_differences
                        ]
                        norm2 = [
                            math.sqrt(sum(dz * dz for dz in diff[1]))
                            for diff in projection_differences
                        ]
                        norminf = [
                            # cannot convert list expr to generator due to
                            #   length check in max_
                            max_([abs(dz) for dz in diff[1]])
                            for diff in projection_differences
                        ]
                        warn(
                            "\tnew points:",
                            ", ".join(map(repr, points_added))
                        )
                        warn("\t1-norm:", ", ".join(map(repr, norm1)))
                        warn("\tsup-norm:", ", ".join(map(repr, norminf)))
                        warn("\t2-norm:", ", ".join(map(repr, norm2)))

                # end of work load discovery

                # save this for status message later before we delete it
                len_extrapolated_items = len(extrapolated_items)

                # turn tobeexplored_boundary and extrapolated_items into
                #   actual list of to be calculated items
                # also setup convergence criterion boundary_projected_count
                boundary_projected_count = float("nan")
                if not work:
                    # Don't let the user interrupt this in the middle
                    #   otherwise the whole projection logic might have to
                    #   be run again
                    with DelayedKeyboardInterrupt():
                        # only calculate something not "nan" if actually
                        #   possible
                        boundary_projected_count = 0
                        # work has been initialized to "set()" before
                        print "# Finding full work load"

                        for index in tobeexplored_boundary:
                            # only explore from "likely" points if
                            #   projections are state-0 complete
                            if not nodes[index][0] >= min_likelihood and \
                               algorithm_state > 1:
                                # use not a >= b instead of < to catch "nan"s
                                # this check might be unnecessary but better
                                #   safe than sorry
                                continue

                            # projection is only done in state 1, 2
                            if algorithm_state in [1, 2] and \
                               is_boundary(index, nodes, par_ranges):
                                boundary_projected_count += 1

                            # nodes[index][1] is 'maybe_boundary' field
                            # if a point is definitely not a boundary point
                            #   there is no need to calculate its neighbors
                            # this check reduces duplicates since points might
                            #   be dropped from nodes but its neighbors still
                            #   remember that they the forgotten ones have
                            #   been calculated
                            if nodes[index][1]:
                                work.update(neighbors(index, par_ranges))

                        work.update(extrapolated_items)

                        # save work load for possible later resuming
                        with open(base_name + ".work", "w") as f:
                            for item in work:
                                f.write(repr(item))
                                f.write("\n")

                        if algorithm_state == 1:
                            # save all projected points
                            with open(base_name + ".projectedpoints", "w") \
                                    as f:
                                for item in tobeexplored_boundary:
                                    f.write("\t".join(map(
                                        out_repr,
                                        index2point(item, par_ranges)
                                    )))
                                    f.write("\n")

                        if cli_arguments.experimental:
                            with open(base_name + ".subnodes", "w") as f:
                                for item in tobeexplored_boundary:
                                    f.write(repr(item))
                                    f.write("\n")

                                # NOTE that items in work can also already
                                #   be in nodes
                                for item in work:
                                    f.write(repr(item))
                                    f.write("\n")

                        del extrapolated_items
                        gc.collect()

                # only calculate new points
                #   e.g. work loaded from resuming could be partially done
                #   already or some neighbors of some boundary points have
                #   been calculated before
                work_items = [
                    (
                        index2point(index, par_ranges),
                        index
                    ) for index in work if index not in nodes
                ]
                work_items.sort(key=lambda p_idx: p_idx[1])

                # tobeexplored_boundary is only available if this is an
                #   iteration where we actually determined the work load
                #   ourselves
                if not tobeexplored_boundary:
                    warn(
                        "# Calculating in total",
                        len(work_items),
                        "new points"
                    )
                else:
                    warn(
                        "# Calculating in total",
                        len(work_items),
                        "new points",
                        "(neighbors of %i points + %i extrapolated points)" % (
                            len(tobeexplored_boundary),
                            len_extrapolated_items
                        )
                    )

                # get rid of the work load in memory here, but keep the file
                #   until the full calculation has been done
                work = set()
                gc.collect()

                category_count = [0, 0, 0]
                # only save progress info over the last 1000 point bunches
                times = deque(
                    [time.time()],
                    maxlen=1000 * concurrent_processors
                )
                full_count = len(work_items)
                for (data, index) in calculate_items(work_items):
                    # point handling is atomic due to delayed interrupts
                    with DelayedKeyboardInterrupt():
                        times.append(time.time())
                        # with open(base_name + ".speed-debug", "a") as f:
                        #     f.write(repr(times[-1]))
                        #     f.write("\n")

                        dressed_data = dress_data(data) \
                            if not isinstance(data, Exception) else None
                        nodes[index] = dressed_data
                        if dressed_data is None:
                            category_count[0] += 1
                        elif not dressed_data[0] >= min_likelihood:
                            # use "not a >= b" instead of "a < b" to catch
                            #   "nan"s
                            category_count[1] += 1
                        else:
                            category_count[2] += 1

                        if not isinstance(data, Exception):
                            suffix = ".data"
                            # also append likelihood to data values
                            data[2].append(dressed_data[0])
                        else:
                            suffix = ".excluded-data"
                            data = [
                                [ERROR_MARKER],
                                list(data.args[0]),
                                [data.args[1]]
                            ]

                        with open(base_name + suffix, "a") as f:
                            f.write("\t".join(map(out_repr, list_sum(data))))
                            f.write("\n")

                        all_category_count = sum(category_count)

                        rate, rateerror, minrate = times_to_rate_info(
                            times,
                            concurrent_processors
                        )

                        draw_progress(
                            all_category_count, full_count,
                            ("fail, unlikely, fine = %i, %i, %i ~ %.1f%%, "
                             "%.1f%%, %.1f%%\nspeed: %f +- %f / s, "
                             "ETA: %s to %s") %
                            tuple(
                                category_count +
                                [
                                    c / (0.01 * all_category_count)
                                    for c in category_count
                                ] + [
                                    rate, rateerror,
                                    progress_forecast(
                                        all_category_count, full_count, rate
                                    ),
                                    progress_forecast(
                                        all_category_count,
                                        full_count, minrate
                                    )
                                ]
                            )
                        )

                # now we can also get rid of the work load file
                #   truncate file by opening it in write mode
                open(base_name + ".work", "w").close()
                # work_items and tobeexplored_boundary may still be needed
                #   for boundary calculation in the next iteration

                valid_count = it_len(
                    k for k, v in nodes.iteritems() if v is not None
                )
                full_count = len(nodes)

                # Note: %s also works if boundary_projected_count is "nan"
                if algorithm_state < 3:
                    # projection mode => tobeexplored_boundary is set of
                    #   projected points
                    warn(
                        "# Point counts: projected&boundary : valid : full",
                        "= %s : %i : %i = %.2f : %.2f : 1" % (
                            boundary_projected_count, valid_count, full_count,
                            boundary_projected_count / full_count,
                            valid_count / full_count
                        )
                    )
                else:
                    # base mode (no projections) => we have full boundary
                    #   available
                    boundary_count = len(boundary)
                    boundary_max_likelihood = max(
                        nodes[index][0] for index in boundary
                    ) if boundary else float("nan")

                    warn(
                        "# Point counts: boundary : valid : full =",
                        "%s : %i : %i = %.2f : %.2f : 1.00" % (
                            boundary_count, valid_count, full_count,
                            boundary_count / full_count,
                            valid_count / full_count
                        )
                    )
                    warn(
                        "# Maximum likelihood on boundary:",
                        boundary_max_likelihood
                    )

                prev_algorithm_state = algorithm_state
                num_explored = len(work_items)

                if work_items:
                    # these new points could open up some paths in the
                    #   projections => reset state
                    algorithm_state = min(available_states)
                    state_one_completeness = 0
                else:
                    if algorithm_state == available_states[-1]:
                        warn(
                            "# State %s exploration yields no new points" %
                            algorithm_state_names[algorithm_state - 1],
                            "- maximal state complete -> scan complete"
                        )
                        break

                    if algorithm_state == 1:
                        warn(
                            "# State one exploration yields no new points%s" %
                            (
                                "(%i times in a row)" % (
                                    state_one_completeness + 1
                                )
                                if state_one_completeness else ""
                            )
                        )
                        state_one_completeness += 1
                    elif algorithm_state == 2:
                        warn(
                            "# State two complete: projections yield no",
                            "new points at all"
                        )
                    else:
                        # this will never be reached as state 3 is global max
                        #   and will be handled by the break in the previous
                        #   conditional block
                        pass

                    algorithm_state = available_states[
                        available_states.index(algorithm_state) + 1
                    ]

                warn(
                    "# Projection algorithm state: %s -> %s" % (
                        algorithm_state_names[prev_algorithm_state - 1],
                        algorithm_state_names[algorithm_state - 1]
                    )
                )

                first_iteration = False

        except Exception as e:
            if cli_arguments.debug and not isinstance(e, FormulaError):
                traceback.print_exc()
            # in case of errors, abandon all stations and let the user fix it
            exit_program("Error: " + str(e))

