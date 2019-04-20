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
Point processor for combination of two or more other point processors.

Usage as script:
    python ProcessorChain.py SCANDEFFILE

Show contents of SCANDEFFILE as seen by the different point processors defined
therein.
"""

import imp
import ConfigParser
processor_modules = []


def init(config_dir, config, module):
    """Initialize processor chain with given processor paths."""
    global processor_modules
    processor_modules = []
    processors = config.get("ProcessorChain", "point_processors").split(":")
    for i, p in enumerate(processors):
        print "# Initializing point processor %i: %s" % (i, p)
        try:
            submodule = imp.load_source(
                "processor_module_%i" % i,
                module.find_file(p, config_dir)
            )
            submodule.init(
                config_dir,
                get_sub_config(i, config),
                module
            )
            processor_modules.append(submodule)
        except IOError as e:
            module.exit_program(
                "Error: could not open point processors '%s': %s" % (p, e)
            )


def main(template_file, pars, vars):
    """Run given point processors in sequence and concatenate result values."""
    global processor_modules
    result = []
    for pm in processor_modules:
        if pm.main.func_code.co_argcount == 1:
            result.extend(pm.main(template_file))
        else:
            result.extend(pm.main(template_file, pars, vars))

    return result


def dict_to_config(d, optionxform):
    """Return ConfigParser instance containing specified settings."""
    config = ConfigParser.SafeConfigParser()
    config.optionxform = optionxform
    for section, items in d.iteritems():
        config.add_section(section)
        for k, v in items.iteritems():
            config.set(section, k, v)
    return config


def get_sub_config(idx, all_config):
    """Return merged config of common sections and specific ones."""
    _dict = all_config._dict
    all_options = _dict(
        (section, _dict(all_config.items(section, raw=True)))
        for section in all_config.sections()
    )
    prefix = "ProcessorChain:%i:" % idx
    translation_table = dict()
    common = _dict(
        s for s in all_options.iteritems()
        if not s[0].startswith("ProcessorChain:")
    )
    specific = _dict(
        (s[0][len(prefix):], s[1]) for s in all_options.iteritems()
        if s[0].startswith(prefix)
    )
    for section, options in common.iteritems():
        for opt in options:
            translation_table[section, opt] = (section, opt)
    for section, options in specific.iteritems():
        for opt in options:
            translation_table[section, opt] = (prefix + section, opt)

    # import pprint
    # pprint.pprint(translation_table)
    merged = _dict(common)
    for section, items in specific.iteritems():
        if section in merged:
            merged[section].update(items)
        else:
            merged[section] = items

    # we use an entire separate ConfigParser object and no proxy/wrapper for
    #   config.get for this, so that interpolation works right
    config = dict_to_config(merged, all_config.optionxform)
    config.get = lambda section, option, raw=False, vars=None, \
        getter=config.get: [
            # first `get' from original config to mark it as used
            all_config.get(
                translation_table[section, option][0],
                translation_table[section, option][1],
                raw=raw, vars=vars
            ),
            # but use the actual merged config, since text interpolation might
            #   make problems otherwise
            getter(section, option, raw, vars)
        ][1]
    return config


if __name__ == "__main__":
    import sys
    if "--help" in sys.argv:
        print __doc__
        sys.exit()

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str  # NOTE keep this in sync with main file
    with open(sys.argv[1]) as f:
        config.readfp(f)

    processors = config.get("ProcessorChain", "point_processors").split(":")
    for i in range(len(processors)):
        subconfig = get_sub_config(i, config)
        print "#" * 40
        print "# As seen by processor %i" % i
        print "#" * 40
        subconfig.write(sys.stdout)
