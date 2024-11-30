import argparse
import os

parser = argparse.ArgumentParser(description='Read and resolve DUNE dependencies given in csv file.')
parser.add_argument('--outputbasedir', nargs = 1, help = 'Directory where to execute "git clone" for all dependencies', required = True)

# parse arguments
args = parser.parse_args()
outputbasedir = args.outputbasedir[0]

with open(f'{os.path.dirname(__file__)}/dune_dependencies.csv') as file:
  dependencies = file.readlines()

# omit first line
dependencies = dependencies[1:]

for dependency in dependencies:
  name, repository, branch, commithash = dependency.split(',')
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
