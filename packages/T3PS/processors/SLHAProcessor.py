#!/usr/bin/env python
# coding=utf-8
#
# Copyright 2014 Vinzenz Maurer
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
Processor that handles SLHA output files and returns their information.

Usage as script:
    python SLHAProcessor.py SLHAFILE [SLHAFILE]...

List parsed information extracted from the SLHAFILEs.
"""

# enabled automatic float division for integers, i.e. 1/2 == 0.5 instead of 0
from __future__ import division

from collections import Sequence
import subprocess
import shlex
import os.path
import os
import warnings
import tempfile

# specifies which column should be treated as value in special blocks
# block names have to be upper case
# for FLHA blocks, if multiple values are given, which differ only by scale,
#   only the last one will be available, as the scale is assumed to be part of
#   the value
# definition format:
#   "BLOCKNAME" : [column indices that give values]
non_standard_blocks = {
    #
    # HiggsBounds BLOCKs
    #

    # format: CoupSq NP IP1 IP2 IP3 ...
    "HIGGSBOUNDSINPUTHIGGSCOUPLINGSBOSONS": [0],

    # format: ScalarNormEffCoupSq PseudoSNormEffCoupSq NP IP1 IP2 IP3 ...
    "HIGGSBOUNDSINPUTHIGGSCOUPLINGSFERMIONS": [0, 1],

    #
    # FLHA BLOCKs
    #

    "FOBS": [2, 3],       # format: ParentPDG type value q NDA ID1 ID2 ID3 ...
    "FOBSSM": [2, 3],     # format: ParentPDG type value q NDA ID1 ID2 ID3 ...

    # format: ParentPDG type minusUncert plusUncert q NDA ID1 ID2 ID3 ...
    "FOBSERR": [2, 3, 4],

    # format: PDG_code mass scheme scale
    "FMASS": [1, 3],

    # format: PDG_code1 code2 nb1 nb2 ratio scheme scale
    "FCONSTRATIO": [4, 6],

    # format: ParentPDG number value NDA ID1 ID2 ID3 ...
    "FPARAM": [2],
}

call_program = data_fields_code = None
slha_files = []
formula_eval = find_binary_file = find_file = timelimit = None
timeout = 10


def init(config_dir, config, module):
    """Initialize SLHA processor."""
    global call_program, slha_files, data_fields_code, timeout, timelimit
    global formula_eval, find_binary_file, find_file
    formula_eval = module.formula_eval
    find_binary_file = module.find_binary_file
    find_file = module.find_file
    timelimit = module.TimeLimit

    print "# Command:", config.get("SLHAProcessor", "program")
    print "#    relative to", config_dir
    try:
        call_program, used_binaries = compile_process_command(
            config.get("SLHAProcessor", "program"),
            config_dir
        )
    except Exception as e:
        module.exit_program("Error in `SLHAProcessor.program': " + str(e))

    print "# Used binaries:", used_binaries

    # get to be read SLHA files from config
    slha_files = config.get("SLHAProcessor", "slha_files").split(":")

    # compile code used to extract interesting information from the SLHA data
    data_fields_code = "(" + config.get("SLHAProcessor", "slha_data") + ")"
    # "eval" mode requires multi-line code to end in "\n"
    if "\n" in data_fields_code and not data_fields_code.endswith("\n"):
        data_fields_code += "\n"

    print "# Extracting fields:", data_fields_code.replace("\n", " ")
    print "#    from files:", slha_files
    try:
        compile(
            data_fields_code,
            data_fields_code,
            "eval",
            flags=division.compiler_flag
        )
    except Exception as e:
        module.exit_program("Error in `SLHAProcessor.slha_data': " + str(e))

    if config.has_option("SLHAProcessor", "timeout"):
        timeout = config.getint("SLHAProcessor", "timeout")
    else:
        timeout = 10
    print "# Timeout:", repr(timeout)


def is_listlike(x):
    """Return True if x is a sequence but not a string."""
    return isinstance(x, Sequence) and not isinstance(x, basestring)


def main(template, pars, vars):
    """Run configured program list and analyze specified SLHA files."""
    global call_program, slha_files, data_fields, timeout, timelimit
    with tempfile.TemporaryFile() as error_sink:
        try:
            with timelimit(timeout):
                if template and os.path.isfile(template):
                    call_program(template, error_sink)
                else:
                    call_program(os.devnull, error_sink)

                slha_data = [read_slha(f) for f in slha_files]

                data = formula_eval(
                    data_fields_code,
                    pars=pars,
                    vars=vars,
                    slha=slha_data
                )

            if is_listlike(data):
                return list(data)
            else:
                return [data]
        except subprocess.CalledProcessError, e:
            error_sink.seek(0)
            e.output = error_sink.read()
            raise e


class SLHABlocks(dict):

    """Dictionary for SLHA file BLOCKs.

    Provides case-insensitive access to BLOCKs, even with scale Q specified.
    """

    # import the global information, so it appears for "help(slha)" and
    #   document it again
    non_standard_blocks = non_standard_blocks
    """Gives the mapping: BLOCK name -> columns that treated as index."""

    def matrix(self, name, Q=Ellipsis):
        """Return the (real or complex) matrix defined as SLHA BLOCK."""
        re = self[name, Q]
        try:
            im = self["IM" + name, Q]
            is_complex = True
        except KeyError:
            is_complex = False

        indices = re.keys() + (im.keys() if is_complex else [])
        if not all(is_listlike(idx) for idx in indices):
            raise ValueError("Block %s is not a proper matrix" % name)

        indices1 = sorted(set(idx[0] for idx in indices))
        indices2 = sorted(set(idx[1] for idx in indices))

        return [
            [
                re[i, j] + (1j * im[i, j] if is_complex else 0.0)
                for j in indices2
            ] for i in indices1
        ]

    def __getitem__(self, requested_index):
        """Return block matching the requested name and possibly scale.

        requested_index can either be the BLOCK name or a tuple consisting of
        BLOCK name and scale Q.

        If scale Q is specified, the block with the nearest scale to the
        requested one is returned. A warning is printed if both scales differ
        by more than 1% of returned one.

        BLOCK name look-up is case-insensitive.

        Examples:
            slhablock["BLOCK"]
            slhablock["BLOCK", 91.1876]
        """
        # if an index of form (name, Q) is requested, find the block with
        #   that name and the nearest Q
        if not isinstance(requested_index, basestring):
            name, Q = requested_index
            name = name.upper()
            blocks_with_Q = [
                k for k in self.keys() if isinstance(k, tuple) and k[0] == name
            ]
            if not blocks_with_Q:
                if name in self:
                    warnings.warn("No %s blocks with specified Q" % name)
                index = name
            elif Q is Ellipsis:
                if len(blocks_with_Q) > 1:
                    raise ValueError(
                        "Multiple blocks named %s, specify scale Q!" % name
                    )
                index = blocks_with_Q[0]
            else:
                index = sorted(blocks_with_Q, key=lambda k: abs(k[1] - Q))[0]
                if abs(Q - index[1]) > 1e-2 * index[1]:
                    warnings.warn(
                        "Nearest Q=%e farther than 1%% from requested %e" %
                        (
                            index[1], Q
                        )
                    )
            return dict.__getitem__(self, index)
        else:
            try:
                return dict.__getitem__(self, requested_index.upper())
            except KeyError:
                return self[requested_index.upper(), ...]


def read_slha(slha_file):
    """Parse SLHA file as defined in [hep-ph/0311123]."""
    currentblock = None
    currentblockname = None
    currentblocktype = None
    currentQ = None
    type_BLOCK = 0
    type_DECAY = 1
    blocks = SLHABlocks()
    decays = dict()
    lineno = 0

    # SLHA file format:
    # lines begin either with BLOCK, DECAY or a space
    # BLOCKs have form: BLOCK <name> [Q= <scale>]
    # DECAYs have form: DECAY <pdgcode> <width>
    # lines beginning with space are data lines and have multiple whitespace
    #   separated data columns of type int, float or string, with the last
    #   one being the "value" and the rest being the index (except for special
    #   cases defined above)
    with open(slha_file) as f:
        for line_ in f:
            lineno += 1
            # split off comments (if existent), strip trailing whitespace
            # WARNING: never left strip since leading whitespace signifies a
            #   data line
            line = line_.rstrip().split("#")[0]
            # str.split without argument splits along one or more whitespace
            parts = line.split()

            if line.upper().startswith("BLOCK"):  # BLOCK declaration line
                if currentblock is not None:
                    # new block -> save current one
                    if currentblocktype == type_BLOCK:
                        if currentQ is None:
                            blocks[currentblockname] = currentblock
                        else:
                            blocks[currentblockname, currentQ] = currentblock
                    elif currentblocktype == type_DECAY:
                        decays[currentblockname] = currentblock

                currentblock = dict()
                currentblocktype = type_BLOCK

                # parts: 0 -> "BLOCK", 1 -> name[, 2 -> "Q=", 3 -> scale]
                currentblockname = parts[1].upper()
                if len(parts) > 2 and parts[2] == "Q=":
                    currentQ = float(parts[3])
                else:
                    currentQ = None

            elif line.upper().startswith("DECAY"):  # DECAY declaration line
                if currentblock is not None:
                    # new block -> save current one
                    if currentblocktype == type_BLOCK:
                        if currentQ is None:
                            blocks[currentblockname] = currentblock
                        else:
                            blocks[currentblockname, currentQ] = currentblock
                    elif currentblocktype == type_DECAY:
                        decays[currentblockname] = currentblock

                currentblock = dict()
                currentblocktype = type_DECAY
                # parts: 0 -> "DECAY", 1 -> PDG code, 2 -> width
                currentblockname = int(parts[1])
                currentblock["width"] = float(parts[2])

            elif line.startswith(" "):  # data line
                # convert each number to float if it contains an exponent or a
                #   decimal place or to integer if not

                # parts = [
                    # maybe_float(p) if ("." in p or "E" in p)
                    # else maybe_int(p)
                    # for p in line.split()
                # ]

                # parts = [maybe_float(p) for p in line.split()]

                # parts = [
                    # maybe_float(p) if looks_like_float(p)
                    # else maybe_int(p)
                    # for p in line.split()
                # ]

                # this seems to be the fastest way to do it:
                parts = [maybe_number(p) for p in line.split()]

                if currentblocktype == type_BLOCK:
                    # line with n parts -> index = first n-1 parts,
                    #   value = last part
                    if currentblockname in non_standard_blocks:
                        value = [
                            x for i, x in enumerate(parts)
                            if i in non_standard_blocks[currentblockname]
                        ]
                        index = [
                            x for i, x in enumerate(parts)
                            if i not in non_standard_blocks[currentblockname]
                        ]
                        currentblock[maybe_tuple(index)] = maybe_tuple(value)
                    else:
                        currentblock[maybe_tuple(parts[:-1])] = parts[-1]
                elif currentblocktype == type_DECAY:
                    # line with n parts -> value = first part,
                    #   index = last n-1 parts
                    currentblock[tuple(parts[1:])] = parts[0]
                else:
                    raise SyntaxError(
                        "Data line before block start",
                        (slha_file, lineno, 0, line)
                    )
            elif not line:
                pass  # ignore empty lines
            else:
                raise SyntaxError(
                    "Unknown line type",
                    (slha_file, lineno, 0, line)
                )

    # also save the last block
    if currentblock is not None:
        if currentblocktype == type_BLOCK:
            if currentQ is None:
                blocks[currentblockname] = currentblock
            else:
                blocks[currentblockname, currentQ] = currentblock
        elif currentblocktype == type_DECAY:
            decays[currentblockname] = currentblock

    blocks["DECAY"] = decays
    return blocks


def looks_like_float(s):
    """Return whether string looks like it could be a float number."""
    return (s == "INF" or s == "NAN" or "." in s or "E" in s)


def maybe_float(s):
    """Try to convert to float, if not possible return as is."""
    try:
        return float(s)
    except Exception:
        return s


def maybe_int(s):
    """Try to convert to int, it not possible, return as is."""
    try:
        return int(s)
    except Exception:
        return s


def maybe_number(s):
    """Try to convert to int, then float, then return as is."""
    try:
        return int(s)
    except Exception:
        pass
    try:
        return float(s)
    except Exception:
        return s


def maybe_tuple(t):
    """Convert an iterable to a) a tuple or b) its first element.

    If only a sequence of one element is passed, this element is returned.
    """
    if len(t) > 1:
        return tuple(t)
    elif len(t) == 1:
        return t[0]
    else:
        return ()


def compile_process_command(source, folder, file_name="<string>"):
    """Compile basic shell-like command script to Python code.

    Arguments:
        folder      base folder from which binaries are searched
        file_name   script file name used for error messages
    """
    # shlex: shell like language lexer, takes care of quotes and escapes
    lex = shlex.shlex(source, posix=True)
    lex.infile = file_name
    lex.whitespace = " \t\r"
    lex.wordchars += "./$[]{}-!:@,+~?*"

    commands = []
    currentcommand = {"args": [], "in": None, "out": None}
    # these tokens govern functionality
    # "<", ">", "&&" as in normal shell
    # ";" delimits commands in the same way as "&&" due to using check_call
    keytokens = set(["&", "&&", ";", "<", ">", lex.eof])
    try:
        for tok in lex:
            # complete single & to && if possible
            if tok == "&":
                nexttoken = lex.get_token()
                if nexttoken != "&":
                    raise SyntaxError(
                        lex.error_leader() +
                        'Expected "&&" instead of "&"'
                    )
                tok += nexttoken

            # ";", "&&" and "\n" delimit commands
            if tok in [";", "&&", "\n"]:
                if not currentcommand["args"] and \
                   (currentcommand["in"] or currentcommand["out"]):
                    raise SyntaxError(
                        lex.error_leader() +
                        'Empty command'
                    )
                if currentcommand["in"] and \
                   currentcommand["in"] == currentcommand["out"]:
                    raise SyntaxError(
                        lex.error_leader() +
                        'Input and output redirection cannot be the same file'
                    )
                commands.append(currentcommand)
                currentcommand = {"args": [], "in": None, "out": None}
            # "< file" specifies input redirection such that the command reads
            #   from file
            elif tok == "<":
                nexttoken = lex.get_token()
                if nexttoken in keytokens:
                    raise SyntaxError(
                        lex.error_leader() +
                        'Incomplete input redirection'
                    )
                currentcommand["in"] = nexttoken
            # "> file" specifies input redirection such that the command writes
            #   to file
            elif tok == ">":
                nexttoken = lex.get_token()
                if nexttoken in keytokens:
                    raise SyntaxError(
                        lex.error_leader() +
                        'Incomplete output redirection'
                    )
                currentcommand["out"] = nexttoken
            # everything else is a normal argument
            else:
                currentcommand["args"].append(tok)

        # don't forget about the last command
        if currentcommand["args"]:
            commands.append(currentcommand)

    except ValueError, err:
        first_line_of_error = lex.token.splitlines()[0]
        raise SyntaxError(
            lex.error_leader() + str(err) +
            ' following "' + first_line_of_error + '"'
        )

    # if the commands use template_placeholder nowhere, pass the template file
    #   to the STDIN of the first command
    template_placeholder = "{template}"
    cmdline_used = False
    for cmd in commands:
        if template_placeholder in cmd["args"] or \
           template_placeholder == cmd["in"]:
            cmdline_used = True
            break

    if not cmdline_used:
        if commands[0]["in"]:
            raise ValueError(
                "Template file name not used and first command already has " +
                "input redirect set"
            )
        commands[0]["in"] = template_placeholder

    # defining a function that is saved in an extra namespace is faster than
    #   compiling just the code, since python only has to do local look-ups for
    #   variables
    # gives about a factor 2 in execution
    tab = " " * 4
    source = "def call_program(template_path, error_sink):\n"
    if cmdline_used and False:
        source += tab + 'with open(os.devnull, "r+") as devnull:\n'
    else:
        source += tab + 'with open(os.devnull, "r+") as devnull, '
        source += 'open(template_path) as template_file:\n'

    for cmd in commands:
        if not cmd["args"]:
            continue

        binary = find_binary_file(cmd["args"][0], folder)

        if not os.path.isfile(binary):
            raise ValueError("Can't find: %s" % cmd["args"][0])
        cmd["args"][0] = binary

        # build up the with clauses needed for input and output redirects, take
        #   care of proper indentation
        with_clause = []
        if cmd["in"] and cmd["in"] != template_placeholder:
            with_clause.append("open(%r) as infile" % cmd["in"])
        if cmd["out"]:
            with_clause.append("open(%r, 'w') as outfile" % cmd["out"])

        if with_clause:
            with_clause = tab * 2 + "with %s:" % ", ".join(with_clause)
            with_clause += "\n" + tab
        else:
            with_clause = ""

        # place calls to specified commands, use check_call so error codes
        #   directly result in exceptions and abortion of the calculation
        source += with_clause + tab * 2
        source += "subprocess.check_call("
        source += "[%s], stdin=%s, stdout=%s, stderr=error_sink)\n" % (
            ", ".join(
                "template_path" if arg == "{template}" else repr(arg)
                for arg in cmd["args"]
            ),
            (
                "infile"
                if cmd["in"] != template_placeholder
                else "template_file"
            ) if cmd["in"] else "devnull",
            "outfile" if cmd["out"] else "devnull"
        )

    # the generated code only needs subprocess, os and built-ins
    ns = {"subprocess": subprocess, "os": os, "open": open}
    exec source in ns
    return ns["call_program"], [cmd["args"][0] for cmd in commands]

if __name__ == "__main__":
    import sys
    if "--help" in sys.argv:
        print __doc__
        sys.exit()

    import pprint  # pretty print
    # argv[0] is the script filename
    slha = []
    for i, file in enumerate(sys.argv[1:]):
        if len(sys.argv) > 2:
            print "File #%s:" % i
        slha.append(read_slha(file))
        pprint.pprint(slha[-1])

    import code  # code -> interactive console
    sys.displayhook = pprint.pprint  # make everything pretty on the console
    warnings.filterwarnings("always")  # always SHOW warnings
    # do not show code in warnings
    warnings.showwarning = lambda message, category, filename, lineno, \
        file=sys.stderr, line=None, original=warnings.showwarning: \
        original(message, category, filename, lineno, file, "")

    code.interact("""
You can now try out formulas involving the "slha" object.
Note that the "math" module (and others) must be included by hand using an
"import" statement before being usable.
Warning: this is a full Python console environment. Be diligent with what you
are doing.

Type "dir(object)" for the list of attributes and methods of "object" or simply
type "dir()" for a list of all objects and functions in the current context.
Type "help" or "exit" for more information.
""",
                  # believe it or not, this indentation is pep8-compliant
                  local={"slha": slha}
                  )
