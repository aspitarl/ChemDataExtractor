ChemDataExtractor
=================

.. image:: http://img.shields.io/pypi/l/ChemDataExtractor.svg?style=flat-square
    :target: https://github.com/mcs07/ChemDataExtractor/blob/master/LICENSE

ChemDataExtractor is a toolkit for extracting chemical information from the scientific literature.
This fork is in the process of implementing a parser for measured chemical concentrations.


Features
--------

- HTML, XML and PDF document readers
- Chemistry-aware natural language processing pipeline
- Chemical named entity recognition
- Rule-based parsing grammars for property and spectra extraction
- Table parser for extracting tabulated data
- Document processing to resolve data interdependencies

Future Improvements
--------

- Support for table layouts
- Improved documentation

Installation
------------
To install this fork of ChemDataExtractor, simply run in a terminal::

    pip install git+https://github.com/maddenfederico/ChemDataExtractor



Documentation
-------------

Full documentation is available at http://chemdataextractor.org/docs

For easy extraction of measured chemical concentrations as a tab-separated spreadsheet, use extract_concentrations.py in an environment with ChemDataExtractor installed
    
``python [path to extract_concentration.py] [path to full-text article]``
   
    

License
-------

ChemDataExtractor is licensed under the `MIT license`_, a permissive, business-friendly license for open source
software.


.. _`MIT license`: https://github.com/mcs07/ChemDataExtractor/blob/master/LICENSE
