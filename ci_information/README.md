<!--
# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
-->

# CI Information

DUNEuro depends on other DUNE modules. Different versions of DUNEuro target different versions of DUNE modules. The file "dune_dependencies.csv" is supposed to define a setup of DUNE modules that works for the present DUNEuro version. The file is structured as follows.

name,repository,description,commithash

Here, "name" assigns a label for the dependency dependency, "repository" is a link to a git repository, "branch" specifies which branch of the repository is supposed to be cloned, and "commithash" specifies the concrete commit hash of the DUNE module this DUNEuro version targets.

We want add some comments.

1. The purpose of "branch" and "commithash" is to provide one example of a matching dependency version. Generally, some other versions of the corresponding DUNE module will also work (and are even explicitly supported). Please refer to the DUNEuro gitlab page for further information ( https://gitlab.dune-project.org/duneuro/duneuro ).
