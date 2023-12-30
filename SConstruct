# -*- coding: utf-8 -*-
#
# __COPYRIGHT__
#
# Copyright © 2012-2024 Philipp Büttgenbach
#
# This file is part of ulphi, a CAE tool and library for computing
# electromagnetic fields.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from datetime import date

EnsureSConsVersion(4, 5)

ProjectName = 'ulphi'
today = date.today()
ProjectVersion = str(today.year-2010)+'.'+str(today.month)
ProjectDict = {
    '@PROJECT_NAME@': ProjectName,
    '@PROJECT_VERSION@': ProjectVersion
}

SConscript(
    'docs/SConscript',
    variant_dir='build/docs', exports="ProjectDict")
