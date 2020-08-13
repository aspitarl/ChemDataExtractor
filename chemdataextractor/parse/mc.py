# -*- coding: utf-8 -*-
"""
chemdataextractor.parse.mc
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Measured Concentration text parser.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import re

from chemdataextractor.parse.cem import cem
from chemdataextractor.parse.common import slash
from chemdataextractor.utils import first
from ..model import Compound, MeasuredConcentration
from .actions import merge
from .base import BaseParser
from .elements import W, R, Optional

log = logging.getLogger(__name__)


si_prefixes_optional = '[YZEPTGMkhcm\u03BCunpfzyad]?'
values = (R('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'))('value')#Matches floats + scientific notation

molar_unit =  R('[YZEPTGMkhcm\u03BCunpfzyad]?M')
mass_unit =  R('[YZEPTGMkhcm\u03BCunpfzyad]?(mol|g)')#I know that moles are a unit of quantity, don't @ me
volume_unit =  R('[YZEPTGMkhcm\u03BCunpfzya]?L')
units = (molar_unit | mass_unit + slash + Optional(volume_unit))('units').add_action(merge)

mc = (values + units)('mc')
# TODO: figure out how to handle creatinine
stopwords = W('was') | W('of') | W('concentration') | W('concentrations')
left_bracket_parentheses = R('[[(]')
right_bracket_parentheses = R('[)\]]')

mc_phrase_cem_first = (cem + Optional(stopwords).hide() + Optional(left_bracket_parentheses).hide() + mc)('mc_phrase')
mc_phrase_cem_last = (mc + Optional(stopwords).hide() + cem)('mc_phrase')

mc_phrase = mc_phrase_cem_first | mc_phrase_cem_last

class McParser(BaseParser):
    """"""
    root = mc_phrase

    def interpret(self, result, start, end):
        compound = Compound(
            measured_concentrations=[
                MeasuredConcentration(
                    value=first(result.xpath('./mc/value/text()')),
                    units=first(result.xpath('./mc/units/text()'))
                )
            ]
        )
        cem_el = first(result.xpath('./cem'))
        if cem_el is not None:
            compound.names = cem_el.xpath('./name/text()')
            compound.labels = cem_el.xpath('./label/text()')
        yield compound