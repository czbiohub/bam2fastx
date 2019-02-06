#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_bam2fastx
----------------------------------

Tests for `bam2fastx` module.
"""

import pytest


class TestBam2fastx(object):

    def test___init__(self, one_hundred):
        assert 100 == one_hundred
