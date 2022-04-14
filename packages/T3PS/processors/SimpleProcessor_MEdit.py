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
Simple processor that calls one program and returns result values.

Usage as script:
    python SimpleProcessor.py < FILE
or
    COMMAND | python SimpleProcessor.py

Show output of COMMAND or content of FILE with extracted numbers highlighted.
"""

# enabled automatic float division for integers, i.e. 1/2 == 0.5 instead of 0
from __future__ import division

from collections import Sequence
import subprocess
import signal
import os.path
import re
import shlex
import commands
arguments = []
timeout = 10
timelimit = None
formula_eval = None
data_fields_code = None
number_pattern = re.compile(
    # make sure that there are () around the full number part
    r"(?<!\w)([-+]?((\d+\.\d*)|(\.\d+)|(\d+))([eE][-+]?\d+)?)"
)


def init(config_dir, config, module):
    """Initialize simple processor.

    Find the absolute path of the requested callable file and prints out what
    it found.
    """
    print '### --- Inside SimpleProcessor.init --- '

    global arguments, timeout, timelimit, data_fields_code, formula_eval
    timelimit = module.TimeLimit
    formula_eval = module.formula_eval


    arguments = shlex.split(config.get("SimpleProcessor", "program"))
    print "                              "
    print "variable 'arguments' is :"
    print arguments
    print "                       "
    original = arguments[0]
    print "variable 'original' is : "
    print original
    print "                    "
    arguments[0] = module.find_binary_file(arguments[0], config_dir)
    if not os.path.isfile(arguments[0]):
        module.exit_program('Error: no such file found: ' + original)
        
    print "# Running", subprocess.list2cmdline(arguments)

    if config.has_option("SimpleProcessor", "timeout"):
        timeout = config.getint("SimpleProcessor", "timeout")
    else:
        timeout = 10
#-- Commented out by Ciara ----------------------------
    print("# Timeout:", timeout, "second" + ("s" if timeout > 1 else ""))
#------------------------------------------------------

    if config.has_option("SimpleProcessor", "data_values"):
        data_fields_code = config.get("SimpleProcessor", "data_values")
        # check for syntax errors
        compile(data_fields_code, repr(data_fields_code), "eval")
        print("# Data values:", data_fields_code)
        print "                          "
        print("#repr(data_fields_code", repr(data_fields_code))
        print("                                  ")
    else:
        data_fields_code = None
  #      print "# Data values: <all>"


def is_listlike(x):
    """Return True if x is a sequence but not a string."""
    return isinstance(x, Sequence) and not isinstance(x, basestring)


def main(template_file, pars, vars):
    """Run requested command and return list of result values."""
   # print '### --- Inside SimpleProcessor.main --- '
    #print '                                        '
    global arguments, timeout, timelimit, formula_eval, data_fields_code

#-- Below was previously commented out - presumably by David
    print 'arguments', arguments
    print 'timeout', timeout
    print 'timelimit', timelimit
    print 'formula_eval', formula_eval
    print'data_fields_code', data_fields_code
    print 'template_file', template_file
#---------------------------------------------------------
   # with open(os.devnull) as devnull, timelimit(timeout):
   #     output = subprocess.check_output(
   #         arguments + ([template_file] if template_file else []),
   #         stdin=devnull,
   #         stderr=subprocess.STDOUT
   #     )
   # print "                                     "
   # print output
   # print "                                     "
#-- Commented out by Ciara ----------------------------
   #print('command output: ', commands.getstatusoutput('{} {}'.format(arguments[0], template_file)) )
#-------------------------------------------------------------------------------------------------------

    with open(os.devnull) as devnull, timelimit(timeout):
        output = subprocess.check_output(
            arguments + ([template_file] if template_file else []),
            stdin=devnull,
            stderr=subprocess.STDOUT
        )


#-- Below was previous commented out - presumably by David
    print "                               "
    print "I am going to print the contents of the template file now"
    print "                               "
#    with open( template_file, 'r' ) as f:
#        content = f.read()
#        print('content: ', content)
#---------------------------------------------
#    command = arguments + [template_file]
#-- Commented out by Ciara-------------
#    print('command', command)
#    print('template_file 3', template_file)
#-------------------------------------------
#    output = subprocess.Popen( command )
#    output.communicate()

    print('### --- End of SimpleProcessor.main --- ')
	
#    def handler(signum,frame):
#        print "Error Occured",signum
#        raise IOError("Segmentation Fault Occured.")
    
#    try:
#        signal.signal(signal.SIGSEGV,handler)  
#    except IOError as e:
#        print e


#    output = subprocess.call(arguments + [template_file])
#--------------------------------------------------------------------

    # group 0 is the full number match
    # make sure it stays that way when changing number_pattern!
    all_numbers = [float(x[0]) for x in re.findall(number_pattern, output)]

#-- Below was previously commented out - presumably by David
#    print('vars', vars)
#    print('pars', pars)
#    print('all_numbers', all_numbers)
#-----------------------------------------------------------------

    return all_numbers
    


    # - This gets called when running test bin
#-- Below was previously commented out - presumably by David
#    if not data_fields_code:
#        return all_numbers
#----------------------------------------------------------------


if __name__ == "__main__":
    import sys
    if "--help" in sys.argv:
        print __doc__
        sys.exit()

    text = sys.stdin.read()
    print "                              "
    print "About to print out sys.stdin.read() - should be the values the processor recieves from parameterprocesser?"
    print "                                         "
    print text
    print "                                "
    def mark_number(s):
        """Mark numbers in color together with index into list of numbers."""
        mark_number.counter += 1
        return (
            "\x1B[35m" + "[%i/%i]" % (
                mark_number.counter,
                -(mark_number.total - mark_number.counter)
            ) + "\x1B[0m" +
            "\x1B[32m" + s.group(0) + "\x1B[0m"
        )

    mark_number.counter = -1
    mark_number.total = len(re.findall(number_pattern, text))

    print "Legend:",
    print "\x1B[35m" + "[index from front/from back]" + "\x1B[32m" + \
        "matched number" + "\x1B[0m"
    print re.sub(number_pattern, mark_number, text)
