<!--
# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
-->

# CI Information

DUNEuro depends on other DUNE modules, and non-DUNE packages, that need to be installed. Different versions of DUNEuro target different versions of DUNE modules. The file "dune_dependencies.csv" is supposed to define a setup of DUNE modules that works for the present DUNEuro version. The file is structured as follows.

name,repository,description,commithash

Here, "name" assigns a label for the dependency, "repository" is a link to a git repository, "branch" specifies which branch of the repository is supposed to be cloned, and "commithash" specifies the concrete commit hash of the DUNE module this DUNEuro version targets.

We want to add some comments.

1. The purpose of "branch" and "commithash" is to provide one example of a matching dependency version. Generally, some other versions of the corresponding DUNE module will also work (and are even explicitly supported). Please refer to the DUNEuro gitlab page for further information ( https://gitlab.dune-project.org/duneuro/duneuro ).
2. For convenience, DUNEuro offers the option to compile bindings to Python and Matlab. The corresponding repositories are described in "bindings.csv", and can be cloned via "get\_python\_bindings.py" and "get\_matlab\_bindings.py". The structure of "bindings.csv" is identical to the structure of "dune_dependencies.csv".
3. To compile duneuro, you can use the "dunecontrol" utility. This is a shell script contained in dune-common/bin. If your working directory contains the cloned DUNE and DUNEuro repositories, you can start the compilation with
  dune-common/bin/dunecontrol --opts=OPTS\_FILE --builddir=BUILD\_PATH all
Here, OPTS\_FILE is a file containing compilation options. This directory contains an exemplary OPTS\_FILE ("config\_release.opts"). BUILD\_PATH is the directory where the compiled code will be put.
4. To compile the Matlab bindings, you need a Matlab installation on your system. Let PATH denote the root of your Matlab installation. You then need to add the line "-DMatlab\_ROOT\_DIR=PATH" to your OPTS\_FILE.
5. A list of the non-DUNE dependencies required for DUNEuro on Ubuntu 22.04 can be found in "resolve\_non\_dune\_dependencies.sh".
