# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Read and get Python bindings for DUNEuro.')
parser.add_argument('--outputbasedir', nargs = 1, help = 'Directory where to execute "git clone"', required = True)

# parse arguments
args = parser.parse_args()
outputbasedir = args.outputbasedir[0]

with open(f'{os.path.dirname(__file__)}/bindings.csv') as file:
  bindings_description = file.readlines()

# omit first line
bindings_description = bindings_description[1:]

duneuropy_found = False
for description in bindings_description:
  name, repository, branch, commithash = description.split(',')
  if name != "duneuro-py":
    continue
  else:
    duneuropy_found = True
    print('Found duneuro-py in binding file')
    print(f'Cloning dependency: {name}')
    print(f'Cloning from: {repository}')
    print(f'Branch: {branch}')
    print(f'Hash: {commithash}')
    os.system(f'mkdir -p {outputbasedir}/{name}')
    os.system(f'git clone --branch {branch} {repository} {outputbasedir}/{name}')
    os.system(f'cd {outputbasedir}/{name} && git checkout {commithash}')
    print(f'Dependency {name} cloned')
    print()
    print()

if not duneuropy_found:
  print('duneuro-py not found in binding file, Python bindings not cloned')
