#! /usr/bin/python3 -t
# -*- coding: utf-8 -*-
#
# Copyright © 2012-2016 Philipp Büttgenbach
#
# This file is part of ulphi, a CAE tool and library for
# computing electromagnetic fields.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.
#

import sys, re
from argparse import ArgumentParser
from os.path import basename

scan_regexs = {
    'include_tex': re.compile(r"\\include\s*\{(\S+)\}"),
    'input_tex': re.compile(r"\\input\s*\{(\S+\.tex)\}"),
    'input_pdf_tex': re.compile(r"\\input\s*\{(\S+\.pdf_tex)\}"),
    'addplot_csv': re.compile(r"\\addplot[^;]+\{(\S+\.csv)\};"),
    'includegraphics': re.compile(r"\\includegraphics[^\{]*\{([^\}]+)\}")
}

results = dict([ (k, set()) for k in scan_regexs.keys() ])

def scan_file(texfilename, quiet=False):
    """Scan file (recursively) for included files.
    """
    
    try:
        texfile = open(texfilename)
    except OSError:
        if not quiet:
            print(basename(sys.argv[0]) + ": Warning: Skipping " + texfilename,
                  file=sys.stderr)
        return

    coding_declaration = re.search(
        r"-\*-.*coding:\s+(\S+).*-\*-", str.join(" ", texfile.readlines(5)))
    if coding_declaration:
        texfile.close()
        try:
            texfile = open(texfilename,
                           encoding=coding_declaration.group(1))
        except LookupError:
            texfile = open(texfilename)

    texfiletext = str.join(" ", texfile)

    local_results = dict([ (k, set()) for k in scan_regexs.keys() ])

    for k in scan_regexs.keys():
        local_results[k].update(
            scan_regexs[k].findall(texfiletext))

    for m in ('include_tex', 'input_tex'):
        for k in local_results[m]:
            # Add .tex if ommitteted
            if not re.search(r"\.[^\.]+$", k):
                local_results[m].remove(k)
                k += ".tex"
                local_results[m].add(k)
            # Also scan that file
            scan_file(k, quiet)

    # Add local results to global results
    for k in results.keys():
        results[k].update(local_results[k])
            

def print_set(s):
    """Sort and print out a set.
    """
    for i in sorted(s):
        print(i, end=";")


def main():
    """Do what must be done.
    """
    option_help_text = {
        'include_tex': "print out \\include{*.tex} files",
        'input_tex': "print out \\input{*.tex} files",
        'input_pdf_tex': "print out \\input{*.pdf_tex} files",
        'addplot_csv': "print out \\addplot table{*.csv}; files",
        'includegraphics': "print out \\includegraphics{*} files"
    }
    
    o_pars = ArgumentParser(
        description="Scan a tex file (recursively) for other included files.")

    o_pars.add_argument(
        '-q', '--quiet', action="store_true", default=False,
        help = "Don't print warnings")
    
    for k in scan_regexs.keys():
        h = ""
        if k in option_help_text:
            h = option_help_text[k]
        o_pars.add_argument(
            "--" + k, action="store_true", default=False, help = h)
        
    o_pars.add_argument(
        "texfilename", type=str, help="(top level) tex file")

    args = o_pars.parse_args()

    scan_file(args.texfilename, args.quiet)

    for k in results.keys():
        if eval('args.' + k):
            print_set(results[k])
            #print()

if "__main__" == __name__ :
    main()
