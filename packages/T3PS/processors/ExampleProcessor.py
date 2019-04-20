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
"""Point processor that showcases the general structure of one."""

# import the math module for most mathematical functions
import math
twopi = 2 * math.pi

# import the cmath module for handling of complex numbers
import cmath
pihalf = -1j * cmath.log(1j)

# we can also get `read_slha' from the SLHAProcessor(.py), which is located
#   in the same directory as this file
import imp
import os.path
read_slha = imp.load_source(
    "SLHAProcessor",
    os.path.join(os.path.dirname(__file__), "SLHAProcessor.py")
).read_slha


def init(config_dir, config, module):
    """Initialize point processor module.

    Arguments:
        config_dir (string): Directory in which the scan definition file is
                             stored.
        config (ConfigParser): Object giving access to definition directives.
        module (Python module): Module containing all functions and classes of
                                the scanner program

    Returns:
        nothing

    Useful functions and objects:
        config.has_section("section")
        config.has_option("section", "directive")
            Check if scan definition contains a particular section or directive
        config.get("section", "directive")
        config.getint("section", "directive")
        config.getfloat("section", "directive")
            Get directive value as string, integer or floating point number
        config.getboolean("section", "directive")
            Get directive value as boolean value
            This function accepts "1", "yes", "true", and "on" as True and
            "0", "no", "false", and "off" as False

        module.formula_eval("python code", pars, vars, data)
            Evaluate formula as Python code.
        module.exit_program("message")
            Exit main program and write message to STDERR.
        module.find_file("relative path")
            Find file from a relative path with the same search rules as for
            the `@include' statement.
        module.find_binary_file("relative path")
            Same as find_file but also searches the PATH environment variable.
        module.TimeLimit(seconds)
            Generate context that has a maximum running time
            Usage:
            with TimeLimit(1):
                <do something that may take longer than 1 second>
    """
    return


def main(template_path, pars, vars):
    """Do or delegate the actual calculation.

    Arguments:
        template_path (string): path to substituted template file that should
                                be processed
        pars (list): list of parameter values
        vars (list): list of variable values
    Returns:
        sequence of result values (float)

    This function is called with the current working directory being the
    temporary directory containing the substituted template file.
    """
    return []
