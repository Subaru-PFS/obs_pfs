#!/usr/bin/env python
from argparse import ArgumentParser
from pfs.utils.headerFixes import HeaderFixDatabase

parser = ArgumentParser(description="Build header correction files")
parser.add_argument("path", help="Path to which to write correction files")
args = parser.parse_args()
HeaderFixDatabase().write(args.path)
