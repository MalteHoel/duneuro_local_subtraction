# SPDX-FileCopyrightText: Copyright © duneuro contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-duneuro-exception OR LGPL-3.0-or-later
# This script was build for the DUNEuro dependencies on Ubuntu 22.04
apt update
apt upgrade -y
apt install -y \
  build-essential \
  cmake \
  git \
  libeigen3-dev \
  libpython3-dev \
  libsuitesparse-dev \
  libtbb-dev \
  libsuperlu-dev
