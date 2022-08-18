# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'WAM2layers'
copyright = '2022, TU Delft, Wageningen University & Research, Netherlands eScience Center'
author = 'Ruud van der Ent, Imme Benedict, Peter Kalverla, Chris Weijenborg'
release = 'v3.0.0-beta'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx_rtd_theme", "myst_nb"]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = '_static/logo.png'


# notebook options
nb_execution_mode = "off"

html_theme_options = {
    "includehidden": True,
    'collapse_navigation': False,
}

# This modifies the "edit on github" link at the top
html_context = {
  'display_github': True,
  'github_user': 'WAM2layers',
  'github_repo': 'WAM2layers',
  'github_version': 'master'
}
