# -*- coding: utf-8 -*-
"""
test_parse_mc
~~~~~~~~~~~~~

Test measured concentration parser

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import unittest

from lxml import etree

from chemdataextractor.doc.text import Sentence
from chemdataextractor.parse.mc import mc_phrase

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


class TestParseMc(unittest.TestCase):
    maxDiff = None

    def do_parse(self, input, expected):
        s = Sentence(input)
        results = []
        for result, *_ in mc_phrase.scan(s.tagged_tokens):
            results.append(etree.tostring(result, encoding='unicode'))
            log.debug(etree.tostring(result, pretty_print=True, encoding='unicode'))
        self.assertEqual(expected, results)

    # Test: "location + cem + infix words + units + normal range"
    def test_mc1(self):
        s = 'The serum homocysteine level was 256.91 μmol/L (normal range, 4.60–12.44)'
        expected = ['<mc_phrase><medium>serum</medium><cem><name>homocysteine</name></cem><mc><value>256.91</value'
                    '><units>μmol/L</units></mc><normal_range><value>4.60–12.44</value></normal_range></mc_phrase>']
        self.do_parse(s, expected)

    # Test: mc_phrase_list
    def test_mc2(self):
        s = 'the 12 most abundant organic metabolites found in serum are D-glucose (5 mM), cholesterol (5 mM), ' \
            'urea (4 mM), ATP (3 mM), glyceraldehyde (1.5 mM),'
        # Can't format it more nicely because then the test fails
        expected = ['<mc_phrase><medium>serum</medium><cem><name>D-glucose</name></cem><mc><value>5</value><units>mM</units></mc><cem><name>cholesterol</name></cem><mc><value>5</value><units>mM</units></mc><cem><name>urea</name></cem><mc><value>4</value><units>mM</units></mc><cem><name>ATP</name></cem><mc><value>3</value><units>mM</units></mc><cem><name>glyceraldehyde</name></cem><mc><value>1.5</value><units>mM</units></mc></mc_phrase>']
        self.do_parse(s, expected)

    def test_mc3(self):
        s = 'The fasting plasma methionine level was increased to 134.34 μmol/L (normal range, 0–90 μmol/L).'
        # Can't format it more nicely because then the test fails
        expected = ['<mc_phrase><medium>plasma</medium><cem><name>methionine</name></cem><mc><value>134.34</value'
                    '><units>μmol/L</units></mc><normal_range><value>0–90</value></normal_range></mc_phrase>']
        self.do_parse(s, expected)

    # Test: Chemical short names next to measured concentrations
    def test_mc4(self):
        s = 'The highest median concentrations for individual analytes in each family were for enterolactone (298 ' \
            'μg/L), monoethylphthalate (MEP; 83.2 μg/L), and benzophenone-3 (BP3; 14.7 μg/L) '
        expected = ['<mc_phrase><cem><name>enterolactone</name></cem><mc><value>298</value><units>μg/L</units></mc'
                    '><cem><name>monoethylphthalate</name></cem><mc><cem><name>MEP</name></cem><value>83.2</value'
                    '><units>μg/L</units></mc></mc_phrase>',
                    '<mc_phrase><cem><name>benzophenone-3</name></cem><mc><cem><name>BP3</name></cem><value>14.7'
                    '</value><units>μg/L '
                    '</units></mc></mc_phrase>']
        self.do_parse(s, expected)

    def test_mc5(self):
        s = 'Her laboratory evaluation showed 48 μM MMA and 11.3 μM MA in plasma and 206 mmol/mol Cr MMA and 26.3 ' \
            'mmol/mol Cr MA in urine '
        # Can't format it more nicely because then the test fails
        expected = ['<mc_phrase><medium>plasma</medium><cem><name>MMA</name></cem><mc><value>134.34</value'
                    '><units>μmol/L</units></mc><normal_range><value>0–90</value></normal_range></mc_phrase>']
        self.do_parse(s, expected)

if __name__ == '__main__':
    unittest.main()
