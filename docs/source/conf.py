# Configuration file for the Sphinx documentation builder.

from scanometrics import _version_major, _version_minor, __version__

# -- Project information

project = 'scanometrics'
copyright = '2024-2025, SCAN-Lab'
author = 'SCAN-Lab'

release = '%d.%d' % (_version_major,_version_minor)
version = __version__

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]
autosummary_generate = True  # Turn on sphinx.ext.autosummary, from https://stackoverflow.com/questions/2701998/sphinx-autodoc-is-not-automatic-enough/62613202#62613202
intersphinx_mapping = {
    'rtd': ("https://docs.readthedocs.io/en/stable/", None),
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
